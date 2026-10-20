// Layout-2 encoding from the BGEN specification. SPDX-License-Identifier: MIT
#include "bgen/encode.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
struct Failure {
    size_t variant;
    size_t sample;
    char const *message;
};

constexpr size_t nowhere = SIZE_MAX;

size_t add(size_t a, size_t b, size_t v, size_t s) {
    if (b > SIZE_MAX - a) throw Failure{v, s, "size arithmetic overflow"};
    return a + b;
}

size_t mul(size_t a, size_t b, size_t v, size_t s) {
    if (a && b > SIZE_MAX / a) throw Failure{v, s, "size arithmetic overflow"};
    return a * b;
}

// The result is bounded before multiplication; no overflowing intermediate.
size_t choose(size_t n, size_t r, size_t v, size_t s) {
    r = std::min(r, n - r);
    size_t value = 1;
    for (size_t j = 1; j <= r; ++j) {
        size_t numerator = n - r + j;
        size_t divisor = j;
        size_t g = std::gcd(numerator, divisor);
        numerator /= g;
        divisor /= g;
        value /= divisor; // After cancellation, divisor divides the previous result.
        value = mul(value, numerator, v, s);
    }
    return value;
}

bool floating(bgen_encode_type t) {
    return t == BGEN_ENCODE_FLOAT32 || t == BGEN_ENCODE_FLOAT64;
}

size_t itemsize(bgen_encode_type t, size_t v) {
    switch (t) {
        case BGEN_ENCODE_FLOAT32: case BGEN_ENCODE_INT32: return 4;
        case BGEN_ENCODE_FLOAT64: case BGEN_ENCODE_INT64: return 8;
        case BGEN_ENCODE_INT16: return 2;
        case BGEN_ENCODE_INT8: return 1;
        default: throw Failure{v, nowhere, "unsupported input type"};
    }
}

template<class T>
T load(void const *base, size_t index) {
    T value;
    std::memcpy(&value, static_cast<uint8_t const *>(base) + index * sizeof(T), sizeof(T));
    return value;
}

int64_t allele(bgen_layout2_variant const &x, size_t index) {
    switch (x.type) {
        case BGEN_ENCODE_INT8: return load<int8_t>(x.data, index);
        case BGEN_ENCODE_INT16: return load<int16_t>(x.data, index);
        case BGEN_ENCODE_INT32: return load<int32_t>(x.data, index);
        case BGEN_ENCODE_INT64: return load<int64_t>(x.data, index);
        default: throw Failure{nowhere, nowhere, "hard calls require signed integers"};
    }
}

double probability(bgen_layout2_variant const &x, size_t index) {
    return x.type == BGEN_ENCODE_FLOAT32
        ? static_cast<double>(load<float>(x.data, index)) : load<double>(x.data, index);
}

struct Call {
    size_t p;
    bool missing;
    size_t active;
    size_t stored;
};

Call call(bgen_layout2_variant const &x, size_t v, size_t s) {
    size_t p = x.ploidy ? x.ploidy[s] : 0;
    bool missing = x.missing && x.missing[s] != 0;
    if (x.missing && x.missing[s] > 1)
        throw Failure{v, s, "missing must contain only 0 or 1"};
    if (!floating(x.type)) {
        size_t observed = 0;
        bool fill = false, absent = false;
        for (size_t h = 0; h < x.width; ++h) {
            int64_t a = allele(x, s * x.width + h);
            if (a == -2) { fill = true; continue; }
            if (fill) throw Failure{v, s, "hard-call -2 fill must be trailing"};
            if (a < -1 || a >= x.nalleles)
                throw Failure{v, s, "hard-call allele outside [-1, K)"};
            absent = absent || a == -1;
            ++observed;
        }
        if (x.ploidy && p != observed)
            throw Failure{v, s, "explicit ploidy disagrees with trailing fill"};
        p = observed;
        if (!x.missing) missing = absent;
        if (!missing && absent)
            throw Failure{v, s, "nonmissing hard call contains -1"};
    }
    if (p > 63) throw Failure{v, s, "ploidy must be in 0..63"};
    size_t active = x.phased ? mul(p, x.nalleles, v, s)
                            : choose(p + x.nalleles - 1, p, v, s);
    size_t stored = x.phased ? mul(p, x.nalleles - 1, v, s) : active - 1;
    if (floating(x.type) && active > x.width)
        throw Failure{v, s, "insufficient probability width for K and ploidy"};
    return {p, missing, active, stored};
}

void descriptor(bgen_layout2_variant const &x, size_t v) {
    if (!x.nalleles) throw Failure{v, nowhere, "allele count must be in 1..65535"};
    if (x.phased > 1) throw Failure{v, nowhere, "phased must be 0 or 1"};
    if (x.bits < 1 || x.bits > 32) throw Failure{v, nowhere, "bits must be in 1..32"};
    size_t bytes = mul(mul(x.nsamples, x.width, v, nowhere), itemsize(x.type, v), v, nowhere);
    if (bytes > x.data_bytes || (bytes && !x.data))
        throw Failure{v, nowhere, "data buffer is too short or null"};
    if ((x.ploidy && x.ploidy_bytes < x.nsamples) ||
        (x.missing && x.missing_bytes < x.nsamples))
        throw Failure{v, nowhere, "ploidy or missing buffer is too short"};
    if (floating(x.type) && x.nsamples && (!x.ploidy || !x.missing))
        throw Failure{v, nowhere, "probability encoding requires ploidy and missing"};
}

size_t block_size(bgen_layout2_variant const &x, size_t v) {
    descriptor(x, v);
    size_t values = 0;
    for (size_t s = 0; s < x.nsamples; ++s)
        values = add(values, call(x, v, s).stored, v, s);
    size_t packed = add(mul(values, x.bits, v, nowhere), 7, v, nowhere) / 8;
    size_t bytes = add(add(10, x.nsamples, v, nowhere), packed, v, nowhere);
    if (bytes > UINT32_MAX)
        throw Failure{v, nowhere, "layout-2 block exceeds uint32 length"};
    return bytes;
}

size_t measure(bgen_layout2_variant const *xs, size_t n, size_t limit) {
    if (n && !xs) throw Failure{nowhere, nowhere, "null variants"};
    if (n > SIZE_MAX / sizeof(*xs)) throw Failure{nowhere, nowhere, "too many variants"};
    size_t total = 0;
    for (size_t v = 0; v < n; ++v) {
        total = add(total, block_size(xs[v], v), v, nowhere);
        if (total > limit) throw Failure{v, nowhere, "batch exceeds max_output_bytes/capacity"};
    }
    return total;
}

struct Packer {
    uint8_t *output;
    uint64_t buffer = 0;
    unsigned available = 0;
    unsigned bits;

    void put(uint64_t value) {
        buffer |= value << available;
        available += bits;
        while (available >= 8) {
            *output++ = static_cast<uint8_t>(buffer);
            buffer >>= 8;
            available -= 8;
        }
    }
    void finish() { if (available) *output++ = static_cast<uint8_t>(buffer); }
};

struct Quantizer {
    std::vector<uint64_t> integers;
    std::vector<double> fractions;
    std::vector<size_t> order;

    void write(bgen_layout2_variant const &x, size_t start, size_t count,
               size_t v, size_t s, uint64_t scale, Packer &packer) {
        double sum = 0;
        for (size_t j = 0; j < count; ++j) {
            double p = probability(x, start + j);
            if (!std::isfinite(p) || p < 0 || p > 1)
                throw Failure{v, s, "probabilities must be finite and in [0, 1]"};
            sum += p;
        }
        if (std::abs(sum - 1) > 1e-6)
            throw Failure{v, s, "probabilities must sum to 1 (absolute tolerance 1e-6)"};
        integers.resize(count);
        fractions.resize(count);
        order.resize(count);
        uint64_t total = 0;
        for (size_t j = 0; j < count; ++j) {
            double scaled = probability(x, start + j) / sum * static_cast<double>(scale);
            integers[j] = static_cast<uint64_t>(std::floor(scaled));
            fractions[j] = scaled - static_cast<double>(integers[j]);
            order[j] = j;
            total += integers[j];
        }
        if (total > scale || scale - total > count)
            throw Failure{v, s, "quantization precision exhausted"};
        size_t remainder = static_cast<size_t>(scale - total);
        // Only membership in the top remainder matters; ties use original index.
        if (remainder) {
            auto compare = [&](size_t a, size_t b) {
                return fractions[a] == fractions[b] ? a < b : fractions[a] > fractions[b];
            };
            if (remainder < count)
                std::nth_element(order.begin(),
                    order.begin() + static_cast<std::vector<size_t>::difference_type>(remainder),
                    order.end(), compare);
            for (size_t j = 0; j < remainder; ++j) ++integers[order[j]];
        }
        for (size_t j = 0; j + 1 < count; ++j) packer.put(integers[j]);
    }
};

void put_le(uint8_t *output, uint64_t value, unsigned bytes) {
    for (unsigned j = 0; j < bytes; ++j) output[j] = static_cast<uint8_t>(value >> (8 * j));
}

void encode(bgen_layout2_variant const &x, size_t v, uint8_t *output, Quantizer &q) {
    put_le(output, x.nsamples, 4);
    put_le(output + 4, x.nalleles, 2);
    uint8_t low = x.nsamples ? 63 : 0, high = 0;
    output[8 + x.nsamples] = x.phased;
    output[9 + x.nsamples] = x.bits;
    Packer packer{output + 10 + x.nsamples, 0, 0, x.bits};
    uint64_t scale = (uint64_t{1} << x.bits) - 1;
    for (size_t s = 0; s < x.nsamples; ++s) {
        Call c = call(x, v, s);
        low = std::min(low, static_cast<uint8_t>(c.p));
        high = std::max(high, static_cast<uint8_t>(c.p));
        output[8 + s] = static_cast<uint8_t>(c.p | (c.missing ? 128 : 0));
        if (c.missing) {
            for (size_t j = 0; j < c.stored; ++j) packer.put(0);
        } else if (floating(x.type)) {
            if (x.phased) {
                for (size_t h = 0; h < c.p; ++h)
                    q.write(x, s * x.width + h * x.nalleles, x.nalleles, v, s, scale, packer);
            } else {
                q.write(x, s * x.width, c.active, v, s, scale, packer);
            }
        } else if (x.phased) {
            for (size_t h = 0; h < c.p; ++h) {
                size_t a = static_cast<size_t>(allele(x, s * x.width + h));
                for (size_t j = 0; j + 1 < x.nalleles; ++j)
                    packer.put(j == a ? scale : 0);
            }
        } else {
            int64_t alleles[63];
            for (size_t h = 0; h < c.p; ++h) alleles[h] = allele(x, s * x.width + h);
            std::sort(alleles, alleles + c.p);
            size_t rank = 0;
            for (size_t h = 0; h < c.p; ++h) {
                size_t n = static_cast<size_t>(alleles[h]) + h;
                if (n >= h + 1) rank = add(rank, choose(n, h + 1, v, s), v, s);
            }
            for (size_t j = 0; j < c.stored; ++j) packer.put(j == rank ? scale : 0);
        }
    }
    packer.finish();
    output[6] = low;
    output[7] = high;
}

void report(bgen_encode_error *error, size_t v, size_t s, char const *message) {
    if (error) {
        error->variant = v;
        error->sample = s;
        std::snprintf(error->message, sizeof(error->message), "%s", message);
    }
}

template<class F>
int guarded(bgen_encode_error *error, F f) {
    report(error, nowhere, nowhere, "");
    try { f(); return 0; }
    catch (Failure const &e) { report(error, e.variant, e.sample, e.message); }
    catch (std::bad_alloc const &) { report(error, nowhere, nowhere, "encoder allocation failed"); }
    catch (std::length_error const &) { report(error, nowhere, nowhere, "encoder allocation too large"); }
    return -1;
}
}

extern "C" int bgen_layout2_size(bgen_layout2_variant const *variants, size_t nvariants,
    size_t max_output_bytes, size_t *size, bgen_encode_error *error) {
    return guarded(error, [&] {
        if (!size) throw Failure{nowhere, nowhere, "null size output"};
        *size = measure(variants, nvariants, max_output_bytes);
    });
}

extern "C" int bgen_layout2_encode(bgen_layout2_variant const *variants, size_t nvariants,
    uint8_t *output, size_t capacity, uint64_t *offsets, size_t offsets_capacity,
    bgen_encode_error *error) {
    return guarded(error, [&] {
        size_t total = measure(variants, nvariants, capacity);
        if (nvariants == SIZE_MAX || offsets_capacity < nvariants + 1 || !offsets)
            throw Failure{nowhere, nowhere, "offsets buffer is too short or null"};
        if (total && !output) throw Failure{nowhere, nowhere, "null payload output"};
        Quantizer q;
        size_t offset = 0;
        offsets[0] = 0;
        for (size_t v = 0; v < nvariants; ++v) {
            size_t bytes = block_size(variants[v], v);
            encode(variants[v], v, output + offset, q);
            offset += bytes;
            offsets[v + 1] = offset;
        }
    });
}

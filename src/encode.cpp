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
#include "zlib-ng.h"

namespace {
struct Failure {
    size_t variant;
    size_t sample;
    char const *message;
    int status = -1;
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

bool fast_hardcall(bgen_layout2_variant const &x) {
    return x.type == BGEN_ENCODE_INT8 && x.nalleles == 2 && x.bits == 8 &&
        (x.width == 1 || x.width == 2) && !x.ploidy && !x.missing;
}

size_t fast_size(bgen_layout2_variant const &x, size_t v) {
    auto *g = static_cast<int8_t const *>(x.data);
    uint64_t values = 0;
    for (size_t s = 0; s < x.nsamples; ++s) {
        int a = g[s * x.width];
        int b = x.width == 2 ? g[s * 2 + 1] : -2;
        if (a < -2 || a > 1 || b < -2 || b > 1)
            throw Failure{v, s, "hard-call allele outside [-1, K)"};
        if (a == -2 && b != -2)
            throw Failure{v, s, "hard-call -2 fill must be trailing"};
        values += static_cast<unsigned>(a != -2) + static_cast<unsigned>(b != -2);
    }
    uint64_t size = 10 + uint64_t{x.nsamples} + values;
    if (size > UINT32_MAX) throw Failure{v, nowhere, "layout-2 block exceeds uint32 length"};
    return static_cast<size_t>(size);
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
    if (fast_hardcall(x)) return fast_size(x, v);
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

// Keep model/width decisions outside the hot loop shared by both public paths.
template<size_t Width, bool Legacy, unsigned Uniform, bool Phase>
size_t fast_encode_one(int8_t const *g, uint32_t ns, uint8_t *out, size_t v) {
    uint8_t *flags = out + 8;
    uint8_t *prob = out + 10 + ns;
    uint8_t low = 2, high = Legacy ? 1 : 0;
    for (size_t s = 0; s < ns; ++s) {
        int a = g[s * Width], b = Width == 2 ? g[s * Width + 1] : -2;
        if (static_cast<unsigned>(a + 2) > 3 || static_cast<unsigned>(b + 2) > 3)
            throw Failure{v, s, "BGEN encoder: genotype value out of range; biallelic input "
                "expects values in {-2, -1, 0, 1}", BGEN_ENCODE_INVALID_ALLELE};
        if (a == -2) {
            if constexpr (Legacy) {
                throw Failure{v, s, "BGEN encoder: -2 in genotype slot 0 (zero-ploidy / unused sample) "
                    "is not representable in this legacy interface", BGEN_ENCODE_INVALID_PLOIDY};
            } else {
                if (b != -2) throw Failure{v, s, "hard-call -2 fill must be trailing"};
                flags[s] = 0;
                low = 0;
                continue;
            }
        }
        if constexpr (Uniform != 0) {
            if ((b == -2) != (Uniform == 1))
                throw Failure{v, s, "fixed layout-2 encoding requires uniform ploidy; "
                    "this chunk contains mixed ploidy", BGEN_ENCODE_MIXED_PLOIDY};
        }
        if (b == -2) {
            flags[s] = static_cast<uint8_t>(a == -1 ? 129 : 1);
            *prob++ = static_cast<uint8_t>(a == 0 ? 255 : 0);
            if constexpr (Uniform == 0) {
                if (low > 1) low = 1;
                if constexpr (!Legacy) {
                    if (high < 1) high = 1;
                }
            }
        } else {
            if (a < 0 || b < 0) {
                flags[s] = 130;
                prob[0] = 0;
                prob[1] = 0;
            } else {
                flags[s] = 2;
                if constexpr (Phase) {
                    prob[0] = static_cast<uint8_t>(a == 0 ? 255 : 0);
                    prob[1] = static_cast<uint8_t>(b == 0 ? 255 : 0);
                } else {
                    prob[0] = static_cast<uint8_t>(a == 0 && b == 0 ? 255 : 0);
                    prob[1] = static_cast<uint8_t>(a != b ? 255 : 0);
                }
            }
            prob += 2;
            if constexpr (Uniform == 0) high = 2;
        }
    }
    put_le(out, ns, 4);
    put_le(out + 4, 2, 2);
    if (!ns) {
        low = high = Legacy ? 2 : 0;
    } else if constexpr (Uniform != 0) {
        low = high = Uniform;
    }
    out[6] = low;
    out[7] = high;
    out[8 + ns] = Phase;
    out[9 + ns] = 8;
    return static_cast<size_t>(prob - out);
}

template<size_t Width, bool Legacy, unsigned Uniform>
size_t fast_model(int8_t const *g, uint32_t ns, uint8_t phase, uint8_t *out, size_t v) {
    return phase ? fast_encode_one<Width, Legacy, Uniform, true>(g, ns, out, v)
                 : fast_encode_one<Width, Legacy, Uniform, false>(g, ns, out, v);
}

size_t fast_encode(int8_t const *g, size_t width, uint32_t ns, uint8_t phase,
                   uint8_t *out, size_t v, bool legacy, unsigned uniform = 0) {
    if (legacy) {
        if (uniform == 1) return fast_model<2, true, 1>(g, ns, phase, out, v);
        if (uniform == 2) return fast_model<2, true, 2>(g, ns, phase, out, v);
        return fast_model<2, true, 0>(g, ns, phase, out, v);
    }
    return width == 1 ? fast_model<1, false, 0>(g, ns, phase, out, v)
                      : fast_model<2, false, 0>(g, ns, phase, out, v);
}

void encode(bgen_layout2_variant const &x, size_t v, uint8_t *output, Quantizer &q) {
    if (fast_hardcall(x)) {
        (void)fast_encode(static_cast<int8_t const *>(x.data), x.width, x.nsamples,
                         x.phased, output, v, false);
        return;
    }
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
    catch (Failure const &e) {
        report(error, e.variant, e.sample, e.message);
        return e.status;
    }
    catch (std::bad_alloc const &) { report(error, nowhere, nowhere, "encoder allocation failed"); }
    catch (std::length_error const &) { report(error, nowhere, nowhere, "encoder allocation too large"); }
    return -1;
}

size_t genotype_size(uint32_t ns, unsigned p) {
    size_t size = add(10, mul(ns, p + 1, nowhere, nowhere), nowhere, nowhere);
    if (size > UINT32_MAX)
        throw Failure{nowhere, nowhere, "layout-2 block exceeds uint32 length"};
    return size;
}

size_t stored_size(size_t bytes) {
    size_t blocks = bytes ? (bytes - 1) / 65535 + 1 : 1;
    return add(add(bytes, mul(blocks, 5, nowhere, nowhere), nowhere, nowhere),
               6, nowhere, nowhere);
}

size_t fixed_size(uint32_t ns, unsigned p, size_t strings) {
    if (p != 1 && p != 2)
        throw Failure{nowhere, nowhere, "uniform_ploidy must be 1 or 2"};
    size_t compressed = stored_size(genotype_size(ns, p));
    if (compressed > UINT32_MAX - 4)
        throw Failure{nowhere, nowhere, "compressed layout-2 framing exceeds uint32 length"};
    return add(add(28, strings, nowhere, nowhere), compressed, nowhere, nowhere);
}

void buffer(void const *data, size_t bytes, size_t need, char const *message) {
    if (bytes < need || (need && !data)) throw Failure{nowhere, nowhere, message};
}

void hardcall_buffers(int8_t const *g, size_t gb, uint8_t const *phase, size_t pb,
                      size_t nv, uint32_t ns) {
    buffer(g, gb, mul(mul(nv, ns, nowhere, nowhere), 2, nowhere, nowhere),
           "genotypes buffer is too short or null");
    buffer(phase, pb, nv, "phased buffer is too short or null");
}

void stored_zlib(uint8_t *out, uint8_t const *data, size_t bytes) {
    *out++ = 0x78;
    *out++ = 1;
    size_t remaining = bytes, offset = 0;
    do {
        size_t count = std::min(remaining, size_t{65535});
        *out++ = static_cast<uint8_t>(remaining <= 65535);
        put_le(out, count, 2);
        put_le(out + 2, static_cast<uint16_t>(~count), 2);
        out += 4;
        std::memcpy(out, data + offset, count);
        out += count;
        offset += count;
        remaining -= count;
    } while (remaining);
    // zlib-ng supplies a platform-optimized checksum; no private zlib clone.
    uint32_t adler = zng_adler32(1, data, static_cast<uint32_t>(bytes));
    for (unsigned j = 0; j < 4; ++j) out[j] = static_cast<uint8_t>(adler >> (24 - 8 * j));
}

size_t string_length(bgen_encode_strings const &s, size_t v) {
    if (!s.width) return 0;
    auto *data = s.data + v * s.width;
    auto *end = static_cast<uint8_t const *>(std::memchr(data, 0, s.width));
    return end ? static_cast<size_t>(end - data) : s.width;
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

extern "C" int bgen_layout2_hardcalls_legacy(
    int8_t const *genotypes, size_t genotypes_bytes, uint8_t const *phased,
    size_t phased_bytes, size_t nvariants, uint32_t nsamples, uint8_t *output,
    size_t output_bytes, size_t row_stride, uint32_t *lengths, size_t lengths_capacity,
    bgen_encode_error *error) {
    return guarded(error, [&] {
        hardcall_buffers(genotypes, genotypes_bytes, phased, phased_bytes, nvariants, nsamples);
        size_t width = genotype_size(nsamples, 2);
        if (row_stride < width || output_bytes < mul(nvariants, row_stride, nowhere, nowhere))
            throw Failure{nowhere, nowhere, "output buffer is too small", BGEN_ENCODE_BUFFER_TOO_SMALL};
        buffer(output, output_bytes, mul(nvariants, row_stride, nowhere, nowhere),
               "null payload output");
        if (lengths_capacity < nvariants || (nvariants && !lengths))
            throw Failure{nowhere, nowhere, "lengths buffer is too short or null"};
        for (size_t v = 0; v < nvariants; ++v) {
            if (phased[v] > 1) throw Failure{v, nowhere, "phased must be 0 or 1"};
            auto *g = nsamples ? genotypes + v * nsamples * 2 : genotypes;
            lengths[v] = static_cast<uint32_t>(
                fast_encode(g, 2, nsamples, phased[v], output + v * row_stride, v, true));
        }
    });
}

extern "C" int bgen_layout2_fixed_size(uint32_t nsamples, uint8_t uniform_ploidy,
    size_t total_string_length, size_t *size, bgen_encode_error *error) {
    return guarded(error, [&] {
        if (!size) throw Failure{nowhere, nowhere, "null size output"};
        *size = fixed_size(nsamples, uniform_ploidy, total_string_length);
    });
}

extern "C" int bgen_layout2_fixed_encode(bgen_layout2_fixed_batch const *batch,
    uint8_t *output, size_t output_bytes, bgen_encode_error *error) {
    return guarded(error, [&] {
        if (!batch) throw Failure{nowhere, nowhere, "null fixed batch"};
        auto const &b = *batch;
        size_t row = fixed_size(b.nsamples, b.uniform_ploidy, b.total_string_length);
        size_t total = mul(b.nvariants, row, nowhere, nowhere);
        if (output_bytes < total)
            throw Failure{nowhere, nowhere, "output buffer is too small", BGEN_ENCODE_BUFFER_TOO_SMALL};
        buffer(output, output_bytes, total, "null output buffer");
        hardcall_buffers(b.genotypes, b.genotypes_bytes, b.phased, b.phased_bytes,
                         b.nvariants, b.nsamples);
        buffer(b.position, b.position_bytes, mul(b.nvariants, 4, nowhere, nowhere),
               "position buffer is too short or null");
        bgen_encode_strings strings[5] = {b.varid, b.rsid, b.chrom, b.allele1, b.allele2};
        for (auto const &s : strings)
            buffer(s.data, s.bytes, mul(b.nvariants, s.width, nowhere, nowhere),
                   "string buffer is too short or null");
        if (!b.nvariants) return;
        size_t raw_size = genotype_size(b.nsamples, b.uniform_ploidy);
        size_t compressed = stored_size(raw_size);
        std::vector<uint8_t> scratch(genotype_size(b.nsamples, 2));
        for (size_t v = 0; v < b.nvariants; ++v) {
            if (b.phased[v] > 1) throw Failure{v, nowhere, "phased must be 0 or 1"};
            auto *g = b.nsamples ? b.genotypes + v * b.nsamples * 2 : b.genotypes;
            (void)fast_encode(g, 2, b.nsamples, b.phased[v], scratch.data(), v, true,
                              b.uniform_ploidy);
            size_t lengths[5], sum = 0;
            for (size_t f = 0; f < 5; ++f) {
                lengths[f] = string_length(strings[f], v);
                if (lengths[f] > (f < 3 ? UINT16_MAX : UINT32_MAX))
                    throw Failure{v, nowhere, "string exceeds BGEN length prefix"};
                sum = add(sum, lengths[f], v, nowhere);
            }
            if (sum != b.total_string_length)
                throw Failure{v, nowhere, "BGEN encoder: per-variant string byte sum does not match "
                    "total_string_length", BGEN_ENCODE_STRING_LENGTH};
            int32_t position = load<int32_t>(b.position, v);
            if (position < 0) throw Failure{v, nowhere, "position must be non-negative"};
            uint8_t *out = output + v * row;
            for (size_t f = 0; f < 5; ++f) {
                if (f == 3) {
                    put_le(out, static_cast<uint32_t>(position), 4);
                    put_le(out + 4, 2, 2);
                    out += 6;
                }
                unsigned prefix = f < 3 ? 2 : 4;
                put_le(out, lengths[f], prefix);
                out += prefix;
                if (lengths[f]) std::memcpy(out, strings[f].data + v * strings[f].width, lengths[f]);
                out += lengths[f];
            }
            put_le(out, compressed + 4, 4);
            put_le(out + 4, raw_size, 4);
            stored_zlib(out + 8, scratch.data(), raw_size);
        }
    });
}

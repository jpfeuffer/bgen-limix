// Strict NumPy batch interfaces. SPDX-License-Identifier: MIT
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION
#include <numpy/arrayobject.h>
#include <nanobind/nanobind.h>
#include "bgen/encode.h"

#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace nb = nanobind;
using namespace nb::literals;

namespace {
struct BufferTooSmall : std::runtime_error {
    using std::runtime_error::runtime_error;
};

PyArrayObject *array(nb::handle value, char const *name, int ndim, int dtype = -1,
                     bool compatibility = false) {
    if (!PyArray_Check(value.ptr()))
        throw nb::type_error((std::string(name) + " must be a NumPy array").c_str());
    auto *a = reinterpret_cast<PyArrayObject *>(value.ptr());
    if (PyArray_NDIM(a) != ndim)
        throw nb::value_error((std::string(name) + " has wrong dimension").c_str());
    if (dtype >= 0 && PyArray_TYPE(a) != dtype) {
        if (compatibility) throw nb::value_error(("Wrong dtype for " + std::string(name)).c_str());
        throw nb::type_error((std::string(name) + " has incorrect dtype").c_str());
    }
    if (!PyArray_IS_C_CONTIGUOUS(a) || !PyArray_ISALIGNED(a) || !PyArray_ISNOTSWAPPED(a))
        throw nb::value_error((std::string(name) +
            " must be aligned, native-endian and C-contiguous (NPY_ARRAY_IN_ARRAY)").c_str());
    return a;
}

void shape(PyArrayObject *a, npy_intp v, npy_intp s, char const *name) {
    if (PyArray_DIM(a, 0) != v || (PyArray_NDIM(a) > 1 && PyArray_DIM(a, 1) != s))
        throw nb::value_error((std::string(name) + " has incompatible variant/sample dimensions").c_str());
}

bgen_encode_type type(PyArrayObject *a, bool hardcall) {
    int t = PyArray_TYPE(a);
    if (!hardcall) {
        if (t == NPY_FLOAT32) return BGEN_ENCODE_FLOAT32;
        if (t == NPY_FLOAT64) return BGEN_ENCODE_FLOAT64;
        throw nb::type_error("gp/hp must have float32 or float64 dtype");
    }
    if (t == NPY_INT8) return BGEN_ENCODE_INT8;
    if (t == NPY_INT16) return BGEN_ENCODE_INT16;
    if (t == NPY_INT32) return BGEN_ENCODE_INT32;
    if (t == NPY_INT64) return BGEN_ENCODE_INT64;
    throw nb::type_error("genotypes must have signed int8, int16, int32 or int64 dtype");
}

uint8_t const *row(PyArrayObject *a, size_t variant, size_t bytes) {
    if (!a || !PyArray_DATA(a)) return nullptr;
    return static_cast<uint8_t const *>(PyArray_DATA(a)) + variant * bytes;
}

void check_error(int status, bgen_encode_error const &error) {
    if (!status) return;
    std::string message;
    if (error.variant != SIZE_MAX) message += "variant " + std::to_string(error.variant);
    if (error.sample != SIZE_MAX) message += ", sample " + std::to_string(error.sample);
    if (!message.empty()) message += ": ";
    message += error.message;
    if (status == BGEN_ENCODE_BUFFER_TOO_SMALL) throw BufferTooSmall(message);
    if (status == BGEN_ENCODE_MIXED_PLOIDY) {
        PyErr_SetString(PyExc_NotImplementedError, message.c_str());
        throw nb::python_error();
    }
    throw nb::value_error(message.c_str());
}

nb::tuple run(std::vector<bgen_layout2_variant> const &variants, size_t max_output_bytes) {
    if (variants.size() >= static_cast<size_t>(std::numeric_limits<npy_intp>::max()))
        throw nb::value_error("too many variants for offsets array");
    size_t bytes = 0;
    bgen_encode_error error{};
    int status;
    {
        nb::gil_scoped_release release;
        status = bgen_layout2_size(variants.data(), variants.size(), max_output_bytes, &bytes, &error);
    }
    check_error(status, error);
    if (bytes > static_cast<size_t>(PY_SSIZE_T_MAX))
        throw nb::value_error("payload exceeds Python bytes capacity");
    nb::object payload = nb::steal<nb::object>(
        PyBytes_FromStringAndSize(nullptr, static_cast<Py_ssize_t>(bytes)));
    if (!payload.is_valid()) throw nb::python_error();
    npy_intp dims[1] = {static_cast<npy_intp>(variants.size() + 1)};
    nb::object offsets = nb::steal<nb::object>(PyArray_EMPTY(1, dims, NPY_UINT64, 0));
    if (!offsets.is_valid()) throw nb::python_error();
    auto *output = reinterpret_cast<uint8_t *>(PyBytes_AsString(payload.ptr()));
    auto *bounds = static_cast<uint64_t *>(PyArray_DATA(
        reinterpret_cast<PyArrayObject *>(offsets.ptr())));
    {
        nb::gil_scoped_release release;
        status = bgen_layout2_encode(variants.data(), variants.size(), output, bytes,
                                    bounds, variants.size() + 1, &error);
    }
    check_error(status, error);
    return nb::make_tuple(payload, offsets);
}

nb::tuple batch(nb::object allele_counts, nb::object ploidy, nb::object missing,
                nb::object phased, nb::object gp, nb::object hp, nb::object genotypes,
                int bits, size_t max_output_bytes) {
    if (bits < 1 || bits > 32) throw nb::value_error("bits must be in 1..32");
    bool hardcall = !genotypes.is_none();
    auto *k = array(allele_counts, "allele_counts", 1, NPY_UINT16);
    auto *phase = array(phased, "phased", 1, NPY_BOOL);
    auto *p = ploidy.is_none() ? nullptr : array(ploidy, "ploidy", 2, NPY_UINT8);
    auto *m = missing.is_none() ? nullptr : array(missing, "missing", 2, NPY_BOOL);
    auto *g = hardcall ? array(genotypes, "genotypes", 3) : nullptr;
    auto *unphased = gp.is_none() ? nullptr : array(gp, "gp", 3);
    auto *haplotypes = hp.is_none() ? nullptr : array(hp, "hp", 3);
    if (!hardcall && (!p || !m)) throw nb::value_error("ploidy and missing are required");
    npy_intp nv = PyArray_DIM(k, 0);
    npy_intp ns = PyArray_DIM(hardcall ? g : p, 1);
    if (ns > UINT32_MAX) throw nb::value_error("sample count exceeds uint32");
    shape(phase, nv, ns, "phased");
    if (p) shape(p, nv, ns, "ploidy");
    if (m) shape(m, nv, ns, "missing");
    for (auto *a : {g, unphased, haplotypes}) {
        if (a) {
            shape(a, nv, ns, "genotypes/gp/hp");
            (void)type(a, a == g);
        }
    }
    auto const *ks = static_cast<uint16_t const *>(PyArray_DATA(k));
    auto const *phases = static_cast<uint8_t const *>(PyArray_DATA(phase));
    std::vector<bgen_layout2_variant> variants(static_cast<size_t>(nv));
    for (size_t v = 0; v < variants.size(); ++v) {
        auto *a = hardcall ? g : (phases[v] ? haplotypes : unphased);
        if (!a)
            throw nb::value_error(("variant " + std::to_string(v) +
                (phases[v] ? ": phased variant requires hp" : ": unphased variant requires gp")).c_str());
        // NumPy ignores strides on singleton axes when checking contiguity.
        size_t row_bytes = static_cast<size_t>(PyArray_NBYTES(a)) / variants.size();
        variants[v] = {
            row(a, v, row_bytes), row_bytes, static_cast<size_t>(PyArray_DIM(a, 2)),
            row(p, v, static_cast<size_t>(ns)), p ? static_cast<size_t>(ns) : 0,
            row(m, v, static_cast<size_t>(ns)), m ? static_cast<size_t>(ns) : 0,
            static_cast<uint32_t>(ns), ks[v], phases[v], static_cast<uint8_t>(bits), type(a, hardcall)
        };
    }
    return run(variants, max_output_bytes);
}

nb::tuple legacy(nb::object genotypes, nb::object phased, size_t limit) {
    auto *g = array(genotypes, "genotypes", 3, NPY_INT8, true);
    auto *p = array(phased, "phased", 1, NPY_BOOL, true);
    if (PyArray_DIM(g, 2) != 2)
        throw nb::value_error("BGEN encoder expects (V, S, 2) genotypes");
    npy_intp nv = PyArray_DIM(g, 0), ns = PyArray_DIM(g, 1);
    if (PyArray_DIM(p, 0) != nv)
        throw nb::value_error("phased.shape[0] must equal genotypes.shape[0]");
    if (ns > (UINT32_MAX - 10) / 3)
        throw nb::value_error("genotype block exceeds uint32 length");
    size_t width = 10 + 3 * static_cast<size_t>(ns);
    if (width > static_cast<size_t>(NPY_MAX_INTP) ||
        static_cast<size_t>(nv) > std::min(limit, static_cast<size_t>(NPY_MAX_INTP)) / width)
        throw nb::value_error("batch exceeds max_output_bytes");
    npy_intp dims[2] = {nv, static_cast<npy_intp>(width)};
    nb::object out = nb::steal<nb::object>(PyArray_EMPTY(2, dims, NPY_UINT8, 0));
    nb::object lens = nb::steal<nb::object>(PyArray_EMPTY(1, dims, NPY_UINT32, 0));
    if (!out.is_valid() || !lens.is_valid()) throw nb::python_error();
    auto *a = reinterpret_cast<PyArrayObject *>(out.ptr());
    auto *l = reinterpret_cast<PyArrayObject *>(lens.ptr());
    bgen_encode_error error{};
    int status;
    {
        nb::gil_scoped_release release;
        status = bgen_layout2_hardcalls_legacy(
            static_cast<int8_t const *>(PyArray_DATA(g)), static_cast<size_t>(PyArray_NBYTES(g)),
            static_cast<uint8_t const *>(PyArray_DATA(p)), static_cast<size_t>(PyArray_NBYTES(p)),
            static_cast<size_t>(nv), static_cast<uint32_t>(ns),
            static_cast<uint8_t *>(PyArray_DATA(a)), static_cast<size_t>(PyArray_NBYTES(a)), width,
            static_cast<uint32_t *>(PyArray_DATA(l)), static_cast<size_t>(nv), &error);
    }
    check_error(status, error);
    return nb::make_tuple(out, lens);
}

size_t fixed_width(Py_ssize_t ns, int ploidy, Py_ssize_t strings) {
    if (ns < 0 || strings < 0) throw nb::value_error("size arguments must be non-negative");
    if (static_cast<size_t>(ns) > UINT32_MAX) throw nb::value_error("sample count exceeds uint32");
    if (ploidy != 1 && ploidy != 2) throw nb::value_error("uniform_ploidy must be 1 or 2");
    size_t result = 0;
    bgen_encode_error error{};
    check_error(bgen_layout2_fixed_size(static_cast<uint32_t>(ns), static_cast<uint8_t>(ploidy),
                                      static_cast<size_t>(strings), &result, &error), error);
    return result;
}

void fixed(nb::object varid, nb::object rsid, nb::object chrom, nb::object allele1,
           nb::object allele2, nb::object position, nb::object genotypes,
           nb::object phased, nb::object out, int uniform_ploidy, Py_ssize_t total_string_length) {
    auto *g = array(genotypes, "genotypes", 3, NPY_INT8, true);
    if (PyArray_DIM(g, 2) != 2)
        throw nb::value_error("BGEN encoder expects (V, S, 2) genotypes");
    npy_intp nv = PyArray_DIM(g, 0), ns = PyArray_DIM(g, 1);
    auto *p = array(phased, "phased", 1, NPY_BOOL, true);
    auto *pos = array(position, "position", 1, NPY_INT32, true);
    if (PyArray_DIM(p, 0) != nv || PyArray_DIM(pos, 0) != nv)
        throw nb::value_error("all per-variant inputs must share num_variants axis");
    auto *id = array(varid, "varid", 2, NPY_UINT8, true);
    auto *rs = array(rsid, "rsid", 2, NPY_UINT8, true);
    auto *chr = array(chrom, "chrom", 2, NPY_UINT8, true);
    auto *a1 = array(allele1, "allele1", 2, NPY_UINT8, true);
    auto *a2 = array(allele2, "allele2", 2, NPY_UINT8, true);
    PyArrayObject *fields[5] = {id, rs, chr, a1, a2};
    bgen_layout2_fixed_batch b{};
    bgen_encode_strings *strings[5] = {&b.varid, &b.rsid, &b.chrom, &b.allele1, &b.allele2};
    for (size_t f = 0; f < 5; ++f) {
        if (PyArray_DIM(fields[f], 0) != nv)
            throw nb::value_error("all per-variant inputs must share num_variants axis");
        *strings[f] = {static_cast<uint8_t const *>(PyArray_DATA(fields[f])),
                      static_cast<size_t>(PyArray_NBYTES(fields[f])),
                      static_cast<size_t>(PyArray_DIM(fields[f], 1))};
    }
    if (!PyArray_Check(out.ptr())) throw nb::type_error("out must be a NumPy array");
    auto *raw = reinterpret_cast<PyArrayObject *>(out.ptr());
    int ndim = PyArray_NDIM(raw);
    if (ndim != 1 && ndim != 2) throw nb::value_error("out_buf has wrong dimension (expected 1D or 2D)");
    auto *output = array(out, "out_buf", ndim, NPY_UINT8, true);
    if (!PyArray_ISWRITEABLE(output)) throw nb::value_error("out_buf must be writeable");
    size_t width = fixed_width(ns, uniform_ploidy, total_string_length);
    if (width > static_cast<size_t>(NPY_MAX_INTP) ||
        static_cast<size_t>(nv) > static_cast<size_t>(NPY_MAX_INTP) / width)
        throw nb::value_error("fixed output dimensions overflow");
    if (ndim == 2 && (PyArray_DIM(output, 0) != nv ||
                      PyArray_DIM(output, 1) != static_cast<npy_intp>(width)))
        throw nb::value_error("out_buf has wrong dimension: expected (variants, layout2_fixed_size(...))");
    if (static_cast<size_t>(PyArray_NBYTES(output)) < static_cast<size_t>(nv) * width)
        throw BufferTooSmall("out_buf is too small");
    auto *out_ptr = static_cast<uint8_t *>(PyArray_DATA(output));
    // Reject overlapping views before native writes can corrupt borrowed input.
    auto overlap = [&](PyArrayObject *input) {
        uintptr_t start = reinterpret_cast<uintptr_t>(PyArray_DATA(input));
        uintptr_t end = start + static_cast<size_t>(PyArray_NBYTES(input));
        uintptr_t out_start = reinterpret_cast<uintptr_t>(out_ptr);
        uintptr_t out_end = out_start + static_cast<size_t>(PyArray_NBYTES(output));
        if (end > start && out_end > out_start && start < out_end && out_start < end)
            throw nb::value_error("out must not overlap input buffers");
    };
    for (auto *field : fields) overlap(field);
    for (auto *input : {g, p, pos}) overlap(input);
    b.position = static_cast<int32_t const *>(PyArray_DATA(pos));
    b.position_bytes = static_cast<size_t>(PyArray_NBYTES(pos));
    b.genotypes = static_cast<int8_t const *>(PyArray_DATA(g));
    b.genotypes_bytes = static_cast<size_t>(PyArray_NBYTES(g));
    b.phased = static_cast<uint8_t const *>(PyArray_DATA(p));
    b.phased_bytes = static_cast<size_t>(PyArray_NBYTES(p));
    b.nvariants = static_cast<size_t>(nv);
    b.nsamples = static_cast<uint32_t>(ns);
    b.uniform_ploidy = static_cast<uint8_t>(uniform_ploidy);
    b.total_string_length = static_cast<size_t>(total_string_length);
    int status;
    bgen_encode_error error{};
    size_t output_bytes = static_cast<size_t>(PyArray_NBYTES(output));
    {
        nb::gil_scoped_release release;
        status = bgen_layout2_fixed_encode(&b, out_ptr, output_bytes, &error);
    }
    check_error(status, error);
}
}

void register_encoding(nb::module_ &m) {
    if (_import_array() < 0) throw nb::python_error();
    nb::exception<BufferTooSmall>(m, "BufferTooSmall", PyExc_ValueError);
    m.def("encode_layout2_hardcalls_legacy", &legacy, "genotypes"_a, "phased"_a,
        nb::kw_only(), "max_output_bytes"_a = 268435456,
        "Encode biallelic int8 (V,S,2) hard calls into padded uint8 rows and uint32 lengths.");
    m.def("layout2_fixed_size", &fixed_width, "num_samples"_a,
          "uniform_ploidy"_a, "total_string_length"_a);
    m.def("encode_layout2_fixed", &fixed, "varid"_a, "rsid"_a, "chrom"_a,
          "allele1"_a, "allele2"_a, "position"_a, "genotypes"_a, "phased"_a,
          "out"_a, "uniform_ploidy"_a, "total_string_length"_a,
          "Write fixed-stride biallelic full variant records with stored zlib compression.");
    m.def("encode_layout2", [](nb::object allele_counts, nb::object ploidy,
        nb::object missing, nb::object phased, nb::object gp, nb::object hp,
        int bits, size_t max_output_bytes) {
            return batch(allele_counts, ploidy, missing, phased, gp, hp, nb::none(),
                         bits, max_output_bytes);
        }, nb::kw_only(), "allele_counts"_a, "ploidy"_a, "missing"_a, "phased"_a,
        "gp"_a = nb::none(), "hp"_a = nb::none(), "bits"_a = 8,
        "max_output_bytes"_a = 268435456,
        "Encode numeric GP/HP into uncompressed layout-2 payload bytes and uint64 offsets.");
    m.def("encode_layout2_hardcalls", [](nb::object genotypes, nb::object allele_counts,
        nb::object phased, nb::object ploidy, nb::object missing, int bits, size_t max_output_bytes) {
            return batch(allele_counts, ploidy, missing, phased, nb::none(), nb::none(),
                         genotypes, bits, max_output_bytes);
        }, nb::kw_only(), "genotypes"_a, "allele_counts"_a, "phased"_a,
        "ploidy"_a = nb::none(), "missing"_a = nb::none(), "bits"_a = 8,
        "max_output_bytes"_a = 268435456,
        "Encode signed integer hard calls (-1 missing, -2 trailing fill) into layout-2 payloads.");
}

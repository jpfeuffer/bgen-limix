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
PyArrayObject *array(nb::handle value, char const *name, int ndim, int dtype = -1) {
    if (!PyArray_Check(value.ptr()))
        throw nb::type_error((std::string(name) + " must be a NumPy array").c_str());
    auto *a = reinterpret_cast<PyArrayObject *>(value.ptr());
    if (PyArray_NDIM(a) != ndim)
        throw nb::value_error((std::string(name) + " has incorrect number of dimensions").c_str());
    if (dtype >= 0 && PyArray_TYPE(a) != dtype)
        throw nb::type_error((std::string(name) + " has incorrect dtype").c_str());
    if (!PyArray_IS_C_CONTIGUOUS(a) || !PyArray_ISALIGNED(a) || !PyArray_ISNOTSWAPPED(a))
        throw nb::value_error((std::string(name) + " must be aligned, native-endian and C-contiguous").c_str());
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
}

void register_encoding(nb::module_ &m) {
    if (_import_array() < 0) throw nb::python_error();
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

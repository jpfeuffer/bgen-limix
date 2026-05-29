#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
#include "bgen/bgen.h"
}

namespace nb = nanobind;
using namespace nb::literals;

// Helper: steal a numpy array PyObject* into a nanobind object.
// The array is allocated by numpy's own allocator (PyArray_EMPTY),
// so no intermediate new[] and no capsule needed.
static nb::object make_numpy_double_2d(npy_intp rows, npy_intp cols, double** out_ptr) {
    npy_intp dims[2] = {rows, cols};
    PyObject* arr = PyArray_EMPTY(2, dims, NPY_DOUBLE, 0 /* C-contiguous */);
    if (!arr) throw std::bad_alloc();
    *out_ptr = static_cast<double*>(PyArray_DATA((PyArrayObject*)arr));
    return nb::steal<nb::object>(arr);
}

static nb::object make_numpy_float_2d(npy_intp rows, npy_intp cols, float** out_ptr) {
    npy_intp dims[2] = {rows, cols};
    PyObject* arr = PyArray_EMPTY(2, dims, NPY_FLOAT, 0);
    if (!arr) throw std::bad_alloc();
    *out_ptr = static_cast<float*>(PyArray_DATA((PyArrayObject*)arr));
    return nb::steal<nb::object>(arr);
}

static nb::object make_numpy_uint8_1d(npy_intp n, uint8_t** out_ptr) {
    npy_intp dims[1] = {n};
    PyObject* arr = PyArray_EMPTY(1, dims, NPY_UINT8, 0);
    if (!arr) throw std::bad_alloc();
    *out_ptr = static_cast<uint8_t*>(PyArray_DATA((PyArrayObject*)arr));
    return nb::steal<nb::object>(arr);
}

static nb::object make_numpy_bool_1d(npy_intp n, bool** out_ptr) {
    npy_intp dims[1] = {n};
    PyObject* arr = PyArray_EMPTY(1, dims, NPY_BOOL, 0);
    if (!arr) throw std::bad_alloc();
    *out_ptr = static_cast<bool*>(PyArray_DATA((PyArrayObject*)arr));
    return nb::steal<nb::object>(arr);
}

// --- Helper RAII wrappers ---

struct BgenFileHandle {
    struct bgen_file* ptr;
    BgenFileHandle(const std::string& filepath) {
        ptr = bgen_file_open(filepath.c_str());
        if (!ptr)
            throw std::runtime_error("Failed to open BGEN file: " + filepath);
    }
    ~BgenFileHandle() {
        if (ptr)
            bgen_file_close(ptr);
    }
    BgenFileHandle(const BgenFileHandle&) = delete;
    BgenFileHandle& operator=(const BgenFileHandle&) = delete;
};

struct BgenMetafileHandle {
    struct bgen_metafile* ptr;
    BgenMetafileHandle(const std::string& filepath) {
        ptr = bgen_metafile_open(filepath.c_str());
        if (!ptr)
            throw std::runtime_error("Failed to open metafile: " + filepath);
    }
    ~BgenMetafileHandle() {
        if (ptr)
            bgen_metafile_close(ptr);
    }
    BgenMetafileHandle(const BgenMetafileHandle&) = delete;
    BgenMetafileHandle& operator=(const BgenMetafileHandle&) = delete;
};

// --- Genotype result struct ---
// Note: We return genotype results as nb::dict to avoid ndarray ownership issues
// with nanobind's reference_internal policy on struct members.

// --- Variant info struct ---

struct VariantInfo {
    std::string id;
    std::string rsid;
    std::string chrom;
    uint32_t position;
    uint16_t nalleles;
    std::vector<std::string> allele_ids;
    uint64_t offset;
};

struct PartitionResult {
    uint32_t offset;
    std::vector<VariantInfo> variants;
};

// --- Python-exposed BgenFile class ---

class BgenFile {
public:
    BgenFile(const std::string& filepath) : handle_(filepath), filepath_(filepath) {}

    std::string filepath() const { return filepath_; }

    int32_t nvariants() const { return bgen_file_nvariants(handle_.ptr); }

    int32_t nsamples() const { return bgen_file_nsamples(handle_.ptr); }

    bool contain_samples() const { return bgen_file_contain_samples(handle_.ptr); }

    std::vector<std::string> read_samples() {
        struct bgen_samples* samples = bgen_file_read_samples(handle_.ptr);
        if (!samples)
            throw std::runtime_error("Could not read samples from BGEN file.");

        int32_t n = nsamples();
        std::vector<std::string> result;
        result.reserve(static_cast<size_t>(n));

        for (int32_t i = 0; i < n; ++i) {
            const struct bgen_string* s = bgen_samples_get(samples, static_cast<uint32_t>(i));
            result.emplace_back(bgen_string_data(s), bgen_string_length(s));
        }

        bgen_samples_destroy(samples);
        return result;
    }

    void create_metafile(const std::string& filepath, bool verbose = false) {
        int32_t nv = nvariants();
        // Estimate partition count: sqrt(nvariants) clamped to [1, nvariants]
        uint32_t nparts = static_cast<uint32_t>(std::max(1.0, std::sqrt(static_cast<double>(nv))));
        if (nparts > static_cast<uint32_t>(nv))
            nparts = static_cast<uint32_t>(nv);

        struct bgen_metafile* mf =
            bgen_metafile_create(handle_.ptr, filepath.c_str(), nparts, verbose ? 1 : 0);
        if (!mf)
            throw std::runtime_error("Error creating metafile: " + filepath);
        bgen_metafile_close(mf);
    }

    nb::dict read_genotype(uint64_t offset) {
        struct bgen_genotype* gt = bgen_file_open_genotype(handle_.ptr, offset);
        if (!gt)
            throw std::runtime_error("Could not open genotype at offset.");

        int32_t ns = nsamples();
        unsigned ncombs = bgen_genotype_ncombs(gt);

        // Allocate directly in numpy's memory — no intermediate new[] or capsule.
        double* prob_data;
        uint8_t* ploidy_data;
        bool* missing_data;
        auto prob_arr   = make_numpy_double_2d(ns, ncombs, &prob_data);
        auto ploidy_arr = make_numpy_uint8_1d(ns, &ploidy_data);
        auto missing_arr = make_numpy_bool_1d(ns, &missing_data);

        int err = bgen_genotype_read64(gt, prob_data);
        if (err != 0) {
            bgen_genotype_close(gt);
            throw std::runtime_error("Could not read genotype probabilities.");
        }

        bool phased = bgen_genotype_phased(gt);

        for (int32_t i = 0; i < ns; ++i) {
            ploidy_data[i] = bgen_genotype_ploidy(gt, static_cast<uint32_t>(i));
            missing_data[i] = bgen_genotype_missing(gt, static_cast<uint32_t>(i));
        }

        bgen_genotype_close(gt);

        nb::dict result;
        result["probability"] = std::move(prob_arr);
        result["phased"] = phased;
        result["ploidy"] = std::move(ploidy_arr);
        result["missing"] = std::move(missing_arr);
        return result;
    }

    nb::dict read_genotype32(uint64_t offset) {
        struct bgen_genotype* gt = bgen_file_open_genotype(handle_.ptr, offset);
        if (!gt)
            throw std::runtime_error("Could not open genotype at offset.");

        int32_t ns = nsamples();
        unsigned ncombs = bgen_genotype_ncombs(gt);

        // Allocate directly in numpy's memory.
        float* prob_data;
        uint8_t* ploidy_data;
        bool* missing_data;
        auto prob_arr    = make_numpy_float_2d(ns, ncombs, &prob_data);
        auto ploidy_arr  = make_numpy_uint8_1d(ns, &ploidy_data);
        auto missing_arr = make_numpy_bool_1d(ns, &missing_data);

        int err = bgen_genotype_read32(gt, prob_data);
        if (err != 0) {
            bgen_genotype_close(gt);
            throw std::runtime_error("Could not read genotype probabilities.");
        }

        bool phased = bgen_genotype_phased(gt);

        for (int32_t i = 0; i < ns; ++i) {
            ploidy_data[i] = bgen_genotype_ploidy(gt, static_cast<uint32_t>(i));
            missing_data[i] = bgen_genotype_missing(gt, static_cast<uint32_t>(i));
        }

        bgen_genotype_close(gt);

        nb::dict result;
        result["probability"] = std::move(prob_arr);
        result["phased"] = phased;
        result["ploidy"] = std::move(ploidy_arr);
        result["missing"] = std::move(missing_arr);
        return result;
    }

    void close() {
        if (handle_.ptr) {
            bgen_file_close(handle_.ptr);
            handle_.ptr = nullptr;
        }
    }

private:
    BgenFileHandle handle_;
    std::string filepath_;
};

// --- Python-exposed BgenMetafile class ---

class BgenMetafile {
public:
    BgenMetafile(const std::string& filepath) : handle_(filepath), filepath_(filepath) {}

    std::string filepath() const { return filepath_; }

    uint32_t npartitions() const { return bgen_metafile_npartitions(handle_.ptr); }

    uint32_t nvariants() const { return bgen_metafile_nvariants(handle_.ptr); }

    uint32_t partition_size() const {
        uint32_t nv = nvariants();
        uint32_t np = npartitions();
        return (nv + np - 1) / np;  // ceildiv
    }

    PartitionResult read_partition(uint32_t index) {
        const struct bgen_partition* partition =
            bgen_metafile_read_partition(handle_.ptr, index);
        if (!partition)
            throw std::runtime_error("Could not read partition.");

        uint32_t nv = bgen_partition_nvariants(partition);
        std::vector<VariantInfo> variants;
        variants.reserve(nv);

        for (uint32_t i = 0; i < nv; ++i) {
            const struct bgen_variant* v = bgen_partition_get_variant(partition, i);
            VariantInfo vi;
            vi.id = std::string(bgen_string_data(v->id), bgen_string_length(v->id));
            vi.rsid = std::string(bgen_string_data(v->rsid), bgen_string_length(v->rsid));
            vi.chrom = std::string(bgen_string_data(v->chrom), bgen_string_length(v->chrom));
            vi.position = v->position;
            vi.nalleles = v->nalleles;
            vi.offset = v->genotype_offset;

            for (uint16_t j = 0; j < v->nalleles; ++j) {
                vi.allele_ids.emplace_back(
                    bgen_string_data(v->allele_ids[j]),
                    bgen_string_length(v->allele_ids[j]));
            }

            variants.push_back(std::move(vi));
        }

        bgen_partition_destroy(partition);

        uint32_t part_offset = partition_size() * index;
        return PartitionResult{part_offset, std::move(variants)};
    }

    void close() {
        if (handle_.ptr) {
            bgen_metafile_close(handle_.ptr);
            handle_.ptr = nullptr;
        }
    }

private:
    BgenMetafileHandle handle_;
    std::string filepath_;
};

// --- Module definition ---

NB_MODULE(_core, m) {
    m.doc() = "C extension for reading BGEN files (nanobind bindings)";
    if (_import_array() < 0)
        throw std::runtime_error("Failed to initialize numpy C API");

    nb::class_<VariantInfo>(m, "VariantInfo")
        .def_ro("id", &VariantInfo::id)
        .def_ro("rsid", &VariantInfo::rsid)
        .def_ro("chrom", &VariantInfo::chrom)
        .def_ro("position", &VariantInfo::position)
        .def_ro("nalleles", &VariantInfo::nalleles)
        .def_ro("allele_ids", &VariantInfo::allele_ids)
        .def_ro("offset", &VariantInfo::offset)
        .def_prop_ro("id_bytes", [](const VariantInfo& v) {
            return nb::bytes(v.id.data(), v.id.size());
        })
        .def_prop_ro("rsid_bytes", [](const VariantInfo& v) {
            return nb::bytes(v.rsid.data(), v.rsid.size());
        })
        .def_prop_ro("chrom_bytes", [](const VariantInfo& v) {
            return nb::bytes(v.chrom.data(), v.chrom.size());
        })
        .def_prop_ro("allele_ids_bytes", [](const VariantInfo& v) {
            // Comma-join allele IDs and return as bytes, matching old cffi format.
            std::string joined;
            for (size_t i = 0; i < v.allele_ids.size(); ++i) {
                if (i > 0) joined += ',';
                joined += v.allele_ids[i];
            }
            return nb::bytes(joined.data(), joined.size());
        });

    nb::class_<PartitionResult>(m, "PartitionResult")
        .def_ro("offset", &PartitionResult::offset)
        .def_ro("variants", &PartitionResult::variants);

    nb::class_<BgenFile>(m, "BgenFile")
        .def(nb::init<const std::string&>(), "filepath"_a)
        .def_prop_ro("filepath", &BgenFile::filepath)
        .def_prop_ro("nvariants", &BgenFile::nvariants)
        .def_prop_ro("nsamples", &BgenFile::nsamples)
        .def_prop_ro("contain_samples", &BgenFile::contain_samples)
        .def("read_samples", &BgenFile::read_samples)
        .def("create_metafile", &BgenFile::create_metafile,
             "filepath"_a, "verbose"_a = false)
        .def("read_genotype", &BgenFile::read_genotype, "offset"_a)
        .def("read_genotype32", &BgenFile::read_genotype32, "offset"_a)
        .def("close", &BgenFile::close)
        .def("__enter__", [](nb::object self) -> nb::object { return self; })
        .def("__exit__", [](BgenFile& self, nb::args) {
            self.close();
        });

    nb::class_<BgenMetafile>(m, "BgenMetafile")
        .def(nb::init<const std::string&>(), "filepath"_a)
        .def_prop_ro("filepath", &BgenMetafile::filepath)
        .def_prop_ro("npartitions", &BgenMetafile::npartitions)
        .def_prop_ro("nvariants", &BgenMetafile::nvariants)
        .def_prop_ro("partition_size", &BgenMetafile::partition_size)
        .def("read_partition", &BgenMetafile::read_partition, "index"_a)
        .def("close", &BgenMetafile::close)
        .def("__enter__", [](nb::object self) -> nb::object { return self; })
        .def("__exit__", [](BgenMetafile& self, nb::args) {
            self.close();
        });
}

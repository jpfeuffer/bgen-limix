#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define NPY_TARGET_VERSION NPY_2_0_API_VERSION
#include <numpy/arrayobject.h>

#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/vector.h>

#include <cmath>
#include <cstdint>
#include <cstring>
#include <algorithm>
#include <limits>
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

static nb::object make_batch_array(std::initializer_list<npy_intp> shape, int type,
                                   size_t item_size, void** data) {
    std::vector<npy_intp> dims(shape);
    size_t count = 1;
    for (npy_intp dimension : dims) {
        if (dimension < 0 || (dimension && count >
            static_cast<size_t>(std::numeric_limits<npy_intp>::max()) /
            static_cast<size_t>(dimension)))
            throw nb::value_error("Batch dimensions overflow.");
        count *= static_cast<size_t>(dimension);
    }
    if (count > static_cast<size_t>(std::numeric_limits<npy_intp>::max()) / item_size)
        throw nb::value_error("Batch byte count overflow.");
    PyObject* array = PyArray_EMPTY(static_cast<int>(dims.size()), dims.data(), type, 0);
    if (!array) throw nb::python_error();
    *data = PyArray_DATA(reinterpret_cast<PyArrayObject*>(array));
    return nb::steal<nb::object>(array);
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

    uint32_t nvariants() const { return bgen_file_nvariants(file()); }

    uint32_t nsamples() const { return bgen_file_nsamples(file()); }

    unsigned layout() const { return bgen_file_layout(file()); }

    bool contain_samples() const { return bgen_file_contain_samples(handle_.ptr); }

    std::vector<std::string> read_samples() {
        struct bgen_samples* samples = bgen_file_read_samples(handle_.ptr);
        if (!samples)
            throw std::runtime_error("Could not read samples from BGEN file.");

        uint32_t n = nsamples();
        std::vector<std::string> result;
        result.reserve(static_cast<size_t>(n));

        for (uint32_t i = 0; i < n; ++i) {
            const struct bgen_string* s = bgen_samples_get(samples, static_cast<uint32_t>(i));
            result.emplace_back(bgen_string_data(s), bgen_string_length(s));
        }

        bgen_samples_destroy(samples);
        return result;
    }

    void create_metafile(const std::string& filepath, bool verbose = false) {
        uint32_t nv = nvariants();
        // Estimate partition count using the same formula as cbgen:
        //   m = max(min(128, nvariants), floor(sqrt(nvariants)))
        //   nparts = nvariants / m  (at least 1)
        uint32_t min_variants = 128;
        uint32_t m = std::max(
            std::min(min_variants, static_cast<uint32_t>(nv)),
            static_cast<uint32_t>(std::floor(std::sqrt(static_cast<double>(nv))))
        );
        uint32_t nparts = (m > 0) ? static_cast<uint32_t>(nv) / m : 1;
        if (nparts < 1) nparts = 1;

        struct bgen_metafile* mf =
            bgen_metafile_create(handle_.ptr, filepath.c_str(), nparts, verbose ? 1 : 0);
        if (!mf)
            throw std::runtime_error("Error creating metafile: " + filepath);
        bgen_metafile_close(mf);
    }

    nb::dict read_genotype(uint64_t offset) {
        struct bgen_genotype* gt = bgen_file_open_genotype(file(), offset);
        if (!gt)
            throw std::runtime_error("Could not open genotype at offset.");

        uint32_t ns = nsamples();
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

        for (uint32_t i = 0; i < ns; ++i) {
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
        struct bgen_genotype* gt = bgen_file_open_genotype(file(), offset);
        if (!gt)
            throw std::runtime_error("Could not open genotype at offset.");

        uint32_t ns = nsamples();
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

        for (uint32_t i = 0; i < ns; ++i) {
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

    /** Return the number of genotype combinations for the variant at *offset*
     *  without allocating or filling any probability array.
     */
    unsigned read_ncombs(uint64_t offset) {
        struct bgen_genotype* gt = bgen_file_open_genotype(file(), offset);
        if (!gt)
            throw std::runtime_error("read_ncombs: could not open genotype at offset.");
        unsigned nc = bgen_genotype_ncombs(gt);
        bgen_genotype_close(gt);
        return nc;
    }

    void close() {
        if (handle_.ptr) {
            bgen_file_close(handle_.ptr);
            handle_.ptr = nullptr;
        }
    }

    /**
     * Return ncombs for each offset (header-only scan, no probability decode).
     * @return 1-D numpy uint32 array of length n_offsets.
     */
    nb::object get_ncombs(const std::vector<uint64_t>& offsets) {
        uint32_t n = batch_size(offsets);
        void* data;
        auto array = make_batch_array({n}, NPY_UINT32, sizeof(uint32_t), &data);
        int err = bgen_file_read_ncombs_batch(file(), offsets.data(), n,
                                             static_cast<uint32_t*>(data));
        if (err)
            throw std::runtime_error("get_ncombs: failed to read ncombs");
        return array;
    }

    nb::dict read_genotype_metadata(const std::vector<uint64_t>& offsets) {
        uint32_t n = batch_size(offsets);
        struct bgen_batch_metadata metadata{};
        auto result = metadata_arrays(n, nsamples(), metadata);
        if (bgen_file_read_genotype_metadata(file(), offsets.data(), n, &metadata))
            throw std::runtime_error("Could not read genotype batch metadata.");
        return result;
    }

    nb::dict read_genotypes_batch_result(const std::vector<uint64_t>& offsets,
                                         int64_t max_ncomb = -1, int precision = 64) {
        uint32_t n = batch_size(offsets);
        uint32_t ns = nsamples();
        struct bgen_batch_metadata metadata{};
        auto result = metadata_arrays(n, ns, metadata);
        result["probability"] = read_batch_probabilities(offsets, max_ncomb, precision, &metadata);
        return result;
    }

    nb::object read_genotypes_batch(const std::vector<uint64_t>& offsets,
                                    int64_t max_ncomb = -1) {
        return read_batch_probabilities(offsets, max_ncomb, 64, nullptr);
    }

private:
    nb::object read_batch_probabilities(const std::vector<uint64_t>& offsets,
                                        int64_t max_ncomb, int precision,
                                        struct bgen_batch_metadata* metadata) {
        uint32_t n = batch_size(offsets);
        uint32_t ns = nsamples();
        if (precision != 32 && precision != 64)
            throw nb::value_error("Precision must be either 64 or 32.");
        if (max_ncomb < -1 || max_ncomb > UINT32_MAX)
            throw nb::value_error("max_ncomb must fit uint32.");
        uint32_t width = max_ncomb == -1 ? 0 : static_cast<uint32_t>(max_ncomb);
        if (max_ncomb == -1 && n) {
            std::vector<uint32_t> counts(n);
            if (bgen_file_read_ncombs_batch(file(), offsets.data(), n, counts.data()))
                throw std::runtime_error("Could not scan genotype batch widths.");
            width = *std::max_element(counts.begin(), counts.end());
        }
        void* out;
        auto probabilities = make_batch_array(
            {n, ns, width}, precision == 32 ? NPY_FLOAT : NPY_DOUBLE,
            precision == 32 ? sizeof(float) : sizeof(double), &out);
        int error = precision == 32 ?
            bgen_file_read_genotypes_batch32(file(), offsets.data(), n, width,
                                             static_cast<float*>(out), metadata) :
            bgen_file_read_genotypes_batch64(file(), offsets.data(), n, width,
                                             static_cast<double*>(out), metadata);
        if (error)
            throw std::runtime_error("Could not read genotype batch (invalid offset or max_ncomb).");
        return probabilities;
    }

    struct bgen_file* file() const {
        if (!handle_.ptr)
            throw std::runtime_error("BGEN file is closed.");
        return handle_.ptr;
    }

    uint32_t batch_size(const std::vector<uint64_t>& offsets) const {
        file();
        if (offsets.size() > UINT32_MAX)
            throw nb::value_error("Too many batch offsets (uint32 limit).");
        return static_cast<uint32_t>(offsets.size());
    }

    static nb::dict metadata_arrays(uint32_t n, uint32_t ns,
                                     struct bgen_batch_metadata& metadata) {
        nb::dict result;
        void* data;
        result["phased"] = make_batch_array({n}, NPY_BOOL, sizeof(bool), &data);
        metadata.phased = static_cast<bool*>(data);
        result["ploidy"] = make_batch_array({n, ns}, NPY_UINT8, sizeof(uint8_t), &data);
        metadata.ploidy = static_cast<uint8_t*>(data);
        result["missing"] = make_batch_array({n, ns}, NPY_BOOL, sizeof(bool), &data);
        metadata.missing = static_cast<bool*>(data);
        result["ncombs"] = make_batch_array({n}, NPY_UINT32, sizeof(uint32_t), &data);
        metadata.ncombs = static_cast<uint32_t*>(data);
        result["probability_count"] = make_batch_array({n, ns}, NPY_UINT32, sizeof(uint32_t), &data);
        metadata.probability_count = static_cast<uint32_t*>(data);
        result["nalleles"] = make_batch_array({n}, NPY_UINT16, sizeof(uint16_t), &data);
        metadata.nalleles = static_cast<uint16_t*>(data);
        result["nbits"] = make_batch_array({n}, NPY_UINT8, sizeof(uint8_t), &data);
        metadata.nbits = static_cast<uint8_t*>(data);
        return result;
    }

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

    /** 1 = all biallelic, 0 = at least one multiallelic, 2 = unknown (v04 metafile) */
    uint8_t all_biallelic() const { return bgen_metafile_all_biallelic(handle_.ptr); }

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
        .def_prop_ro("layout", &BgenFile::layout)
        .def_prop_ro("contain_samples", &BgenFile::contain_samples)
        .def("read_samples", &BgenFile::read_samples)
        .def("create_metafile", &BgenFile::create_metafile,
             "filepath"_a, "verbose"_a = false)
        .def("read_genotype", &BgenFile::read_genotype, "offset"_a)
        .def("read_genotype32", &BgenFile::read_genotype32, "offset"_a)
        .def("read_ncombs", &BgenFile::read_ncombs, "offset"_a,
             "Return the number of genotype combinations for the variant at offset "
             "without materialising the probability array.")
        .def("get_ncombs", &BgenFile::get_ncombs, "offsets"_a,
             "Return a uint32 array of ncombs values for each offset (header-only scan).")
        .def("read_genotype_metadata", &BgenFile::read_genotype_metadata, "offsets"_a)
        .def("read_genotypes_batch_result", &BgenFile::read_genotypes_batch_result,
             "offsets"_a, "max_ncomb"_a = -1, "precision"_a = 64)
        .def("read_genotypes_batch", &BgenFile::read_genotypes_batch,
             "offsets"_a, "max_ncomb"_a = -1,
             "Read probabilities for multiple variants into a padded rectangular array. "
             "Returns float64 array of shape (n_variants, nsamples, max_ncomb). "
             "Unused combination slots are NaN-filled. "
             "Pass max_ncomb=-1 (default) to auto-detect from all offsets.")
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
        .def_prop_ro("all_biallelic", &BgenMetafile::all_biallelic,
             "1=all biallelic, 0=has multiallelic, 2=unknown (v04 metafile)")
        .def("read_partition", &BgenMetafile::read_partition, "index"_a)
        .def("close", &BgenMetafile::close)
        .def("__enter__", [](nb::object self) -> nb::object { return self; })
        .def("__exit__", [](BgenMetafile& self, nb::args) {
            self.close();
        });
}

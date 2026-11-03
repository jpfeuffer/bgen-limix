"""Tests for cbgen Python bindings."""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from cbgen import BgenFile, BgenMetafile, example


@pytest.mark.parametrize("legacy", [False, True])
@pytest.mark.parametrize("cloud_object", [False, True])
def test_s3a_path_normalization(monkeypatch, legacy, cloud_object):
    import cbgen
    from cbgen import _compat

    class CloudPath:
        storage_options = {"key": "path-key", "secret": "path-secret"}

        def __str__(self):
            return "s3a://bucket/example.bgen"

    def native_open(path):
        import os

        assert path == "s3://bucket/example.bgen"
        if cloud_object:
            assert os.environ["AWS_ACCESS_KEY_ID"] == "path-key"
            assert os.environ["AWS_SECRET_ACCESS_KEY"] == "path-secret"
        return SimpleNamespace(filepath=path, close=lambda: None)

    monkeypatch.setattr(_compat if legacy else cbgen, "_BgenFile", native_open)
    factory = _compat.bgen_file if legacy else cbgen.BgenFile
    reader = factory(CloudPath() if cloud_object else "s3a://bucket/example.bgen")
    reader.close()


def test_haplotypes_basic():
    filepath = example.get("haplotypes.bgen")
    with BgenFile(str(filepath)) as bgen:
        assert bgen.nvariants == 4
        assert bgen.nsamples == 4
        assert bgen.contain_samples

        samples = bgen.read_samples()
        assert samples == ["sample_0", "sample_1", "sample_2", "sample_3"]


def test_haplotypes_genotype(tmp_path):
    filepath = example.get("haplotypes.bgen")
    mfpath = tmp_path / "haplotypes.bgen.metafile"

    with BgenFile(str(filepath)) as bgen:
        bgen.create_metafile(str(mfpath))
        with BgenMetafile(str(mfpath)) as mf:
            assert mf.npartitions == 1
            assert mf.nvariants == 4
            assert mf.partition_size == 4

            part = mf.read_partition(0)
            assert len(part.variants) == 4

            gt = bgen.read_genotype(part.variants[0].offset)
            expected = np.array([
                [1.0, 0.0, 1.0, 0.0],
                [0.0, 1.0, 1.0, 0.0],
                [1.0, 0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0, 1.0],
            ])
            np.testing.assert_allclose(gt.probability, expected)
            assert gt.phased is True


def test_metafile_creation(tmp_path):
    filepath = example.get("haplotypes.bgen")
    mfpath = tmp_path / "test.metafile"

    with BgenFile(str(filepath)) as bgen:
        bgen.create_metafile(str(mfpath))
        assert mfpath.exists()

        with BgenMetafile(str(mfpath)) as mf:
            assert mf.nvariants == 4


def test_variant_info(tmp_path):
    filepath = example.get("haplotypes.bgen")
    mfpath = tmp_path / "haplotypes.bgen.metafile"

    with BgenFile(str(filepath)) as bgen:
        bgen.create_metafile(str(mfpath))

    with BgenMetafile(str(mfpath)) as mf:
        part = mf.read_partition(0)
        v = part.variants[0]
        assert v.id == "SNP1"
        assert v.rsid == "RS1"
        assert v.chrom == "1"
        assert v.position == 1
        assert v.nalleles == 2
        assert v.allele_ids == ["A", "G"]


def test_context_manager():
    filepath = example.get("haplotypes.bgen")
    with BgenFile(str(filepath)) as bgen:
        assert bgen.nvariants == 4
    # Should not raise after close


def test_no_samples():
    filepath = example.get("complex.23bits.no.samples.bgen")
    with BgenFile(str(filepath)) as bgen:
        assert not bgen.contain_samples
        with pytest.raises(RuntimeError):
            bgen.read_samples()

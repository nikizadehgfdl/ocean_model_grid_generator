# The checksums in this file were generated on gfdl PAN platform with
# commit ca730b83714630a05c4c8f133b23152a6165909a
import pytest
import subprocess
import hashlib
import os
import glob

my_srcdir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ogg_cmd = os.path.join(my_srcdir, "ocean_grid_generator.py")


def hashfile(f):
    sha256 = hashlib.sha256()
    with open(f, "rb") as file:
        sha256.update(file.read())
    return sha256.hexdigest()


class TestOGG:
    def test_hgrid_res4_0(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res4.0.nc")
        sp = subprocess.run(
            [ogg_cmd, "-f", outfile, "-r", "0.25", "--ensure_nj_even", "--no_changing_meta"]
        )
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "47f5fc42e7f598f1d339c15e7fc88dc1d7b03956898294f4f185fe17b0a0f31d"
        )

    def test_hgrid_res1_0(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res1.0.nc")
        sp = subprocess.run(
            [
                ogg_cmd,
                "-f",
                outfile,
                "-r",
                "1.0",
                "--south_cutoff_row",
                "2",
                "--no_changing_meta",
            ]
        )
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "bb57a86b788cb27c71ad139b5c72892fbe7019ffe774ecce4c0d74de54e678dc"
        )

    def test_hgrid_res0_5(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res0.5.nc")
        sp = subprocess.run([ogg_cmd, "-f", outfile, "-r", "2", "--no_changing_meta"])
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "87b29a240aa1ea9dbcce7d718a13704fa96b098650acbd99e4376884c35c4c83"
        )

    def test_hgrid_res0_5_equenh(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res0.5_equenh.nc")
        sp = subprocess.run(
            [
                ogg_cmd,
                "-f",
                outfile,
                "-r",
                "2.0",
                "--south_cutoff_row",
                "130",
                "--no_changing_meta",
                "--write_subgrid_files",
                "--enhanced_equatorial",
            ]
        )
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "683a42d5a155620456635dac57a2a6789058fd2a496feafb05db5b2c6015e754"
        )

    def test_hgrid_res0_25(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res0.25.nc")
        sp = subprocess.run(
            [
                ogg_cmd,
                "-f",
                outfile,
                "-r",
                "4",
                "--r_dp",
                "0.2",
                "--south_cutoff_row",
                "83",
                "--write_subgrid_files",
                "--no_changing_meta",
            ]
        )
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "12185aed9f1814c6b3c1545158338ee4160da96e2c9f389ffc31e1d3ed8c76ca"  # gfdl-pan106,pan105
            or hashfile(outfile)
            == "2d7840344aa356feb282d352bf21cce8832947b7989534290bd24d30dc561b70"  # github,gfdl-pan202
        )


## test_hgrid_res0_125 might require more memory than available on github Action platform
#  def test_hgrid_res0_125(self, tmpdir):
#    outfile = tmpdir.join('ocean_hgrid_res0.125.nc')
#    sp = subprocess.run([ogg_cmd,
#                         '-f', outfile, '-r', '8', '--r_dp', '0.2',
#                         '--south_cutoff_row', '5', '--match_dy', '--ensure_nj_even',
#                         '--write_subgrid_files', '--no_changing_meta'])
#    assert sp.returncode == 0
#    assert hashfile(outfile) == 'a5a5c859464aecda946a8d43ab238256fe26c5bec546a2d8d533de221276a63f'#gfdl-pan106,pan105
#    #assert hashfile(outfile) == 'f1eb1fec9e47110013713c2fc69169a98423644673d9002749f5ded6c4b0a09b' #gfdl-pan202

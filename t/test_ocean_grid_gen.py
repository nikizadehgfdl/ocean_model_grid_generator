# The checksums in this file were generated on gfdl PAN platform with
# commit ca730b83714630a05c4c8f133b23152a6165909a
import pytest
import subprocess
import hashlib
import os
import glob
import sys

my_srcdir = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
ogg_cmd = os.path.join(my_srcdir, "ocean_grid_generator.py")


def hashfile(f):
    sha256 = hashlib.sha256()
    with open(f, "rb") as file:
        sha256.update(file.read())
    return sha256.hexdigest()

# Note, all calls to subprocess use `sys.executable` to ensure the Python
# executable pytest uses is used.  This is to help isolate if the shebang
# in ocean_grid_generatory.py points to a version that does not have
# thre required python modules.
class TestOGG:
    def test_hgrid_res4_0(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res4.0.nc")
        sp = subprocess.run(
            [
                sys.executable,
                ogg_cmd,
                "-f",
                outfile,
                "-r",
                "0.25",
                "--ensure_nj_even",
                "--no_changing_meta"
            ]
        )
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "0bc8e50b7d4110fcfb148f9f9590cf3c5e6c2d699b9c5519600aea7cc9483381" #PAN
            or hashfile(outfile)
            == "7e04457b4411260341f304a5182ec7f289ddb2724021e2825751f9994e23a6ed" #githubCI
        )

    def test_hgrid_res1_0(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res1.0.nc")
        sp = subprocess.run(
            [
                sys.executable,
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
            == "6fed3c3f57894ea5e635c5412fdcdef90f44576908b10a894f06035c38159f9f" #PAN 
            or hashfile(outfile)
            == "ec24bbb7fd8b076d4d566777d2a924c4c6ac8f45bdf5b55c94eaeaa12ab07cbf" #githubCI
        )

    def test_hgrid_res0_5(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res0.5.nc")
        sp = subprocess.run(
            [
                sys.executable,
                ogg_cmd,
                "-f",
                outfile,
                "-r",
                "2",
                "--no_changing_meta"
            ]
        )
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "cadb254616cee8c0fd8ffc52360de7d5b23ef298397e2873151b4c7e55d999b7" #PAN
            or hashfile(outfile)
            == "b2e8bbf65da5aa1de011119a269f4b3f236bf62ebe14778f1f458b552a1c1a9d" #githubCI
        )

    def test_hgrid_res0_5_equenh(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res0.5_equenh.nc")
        sp = subprocess.run(
            [
                sys.executable,
                ogg_cmd,
                "-f",
                outfile,
                "-r",
                "2.0",
                "--no_changing_meta",
                "--write_subgrid_files",
                "--enhanced_equatorial=4",
            ]
        )
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "26cee1438eb2a5188ef708528b53fca07bf2312b6cfd92f8e76dbf7f6c6ecaa0" #PAN
            or hashfile(outfile)
            == "30a040b8fbdf0683cc50d8ac2be841bb81f3f681e3b0c1a274f841e6fbd47d44" #githubCI
        )

    def test_hgrid_res0_25(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res0.25.nc")
        sp = subprocess.run(
            [
                sys.executable,
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
            == "fd8ccd83cf5fff0d6c4a0a3abf7f5a4f561107ddf4f282f1289a113077ca6738"  # gfdl-pan106,pan105
            or hashfile(outfile)
            == "2d7840344aa356feb282d352bf21cce8832947b7989534290bd24d30dc561b70"  # github,gfdl-pan202
            or hashfile(outfile)
            == "1a0c3ca0e5b71ebdb8c68c63b2134c48a4bc8d9b9c9ba32f7298171a315d9508"  # githubCI
            #== "20315d58d80747559ed069fff521511624fed778ab92df6f9764503b79a15eea"  # githubCI old
        )

    def test_hgrid_res0_25_om5proto(self, tmpdir):
        outfile = tmpdir.join("ocean_hgrid_res0.25_om5proto.nc")
        sp = subprocess.run(
            [
                sys.executable,
                ogg_cmd,
                "-f", outfile, "-r","4",
                "--south_ocean_lower_lat", "-88.57", "--match_dy", "so", "--no_south_cap",
                "--write_subgrid_files","--no_changing_meta",
            ]
        )
        assert sp.returncode == 0
        assert (
            hashfile(outfile)
            == "0dd436ac1de339c447bd53a47f54cdff153a73a31e782bee59b1699209aa6b0d"  # gfdl-pan106,pan105
            or hashfile(outfile)
            == "1a0c3ca0e5b71ebdb8c68c63b2134c48a4bc8d9b9c9ba32f7298171a315d9508"  # githubCI
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

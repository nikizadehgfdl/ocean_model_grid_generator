import pytest
import subprocess
import hashlib
import os
import glob

my_srcdir=os.path.abspath(os.path.join(os.path.dirname(__file__),'..'))
ogg_cmd=os.path.join(my_srcdir, 'ocean_grid_generator.py')

def hashfile(f):
  sha256 = hashlib.sha256()
  with open(f, 'rb') as file:
    sha256.update(file.read())
  return sha256.hexdigest()

class TestOGG():

  def test_hgrid_res4_0(self, tmpdir):
    tic0 = time.perf_counter()
    outfile = tmpdir.join('ocean_hgrid_res4.0.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '0.25',
                         '--even_j', '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile) == 'ccf9b9bcb789e57722b8af867079eb3ad6788da98f5fa855aeba029771e116a8'
    print("res4 time(s): ",time.perf_counter() - tic0)

  def test_hgrid_res1_0(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res1.0.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '1.0',
                         '--south_cutoff_row', '2', '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile) == 'fd15681f5fb8ee81b7c9f7aa4f3a51291161b4c77a6594cae46db18043affc48'

  def test_hgrid_res0_5(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res0.5.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '2',
                         '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile) == '2407058e0e2685b52c5cb6d5ca3c94341cacf6aea52c9e588ca08a13d7bbfd75'
  
  def test_hgrid_res0_5_equenh(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res0.5_equenh.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '2.0',
                         '--south_cutoff_row', '130', '--no_changing_meta',
                         '--write_subgrid_files', '--enhanced_equatorial'])
    assert sp.returncode == 0
    assert hashfile(outfile) == '223250ce5aa1f0bd2e3f976e84583a1b5836a4a8323e64119027fbc165e370a2'

  def test_hgrid_res0_25(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res0.25.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '4', '--rdp', '0.2',
                         '--south_cutoff_row', '83', '--write_subgrid_files', '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile) == '299f16cc397db0b1a595fb3d5d81465479452679a77530b6dc48a74e3c30af9b'

  def test_hgrid_res0_125(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res0.125.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '8', '--rdp', '0.2',
                         '--south_cutoff_row', '5', '--match_dy', '--even_j',
                         '--write_subgrid_files', '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile) == 'fed9912406f2fda95fd7ab81a3b7b7c432c4a6834c9ba0f46afceab2e1ad6264'

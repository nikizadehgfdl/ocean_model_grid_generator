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
    outfile = tmpdir.join('ocean_hgrid_res4.0.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '0.25',
                         '--even_j', '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile), '3918a72255449958108dd68e99809537b0f29147b341e06ce5e4c341c855e245'

  def test_hgrid_res1_0(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res1.0.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '1.0',
                         '--south_cutoff_row', '2', '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile), 'f25591f7b905527c744ed568ca3e211c4be803e69d313acfe0c697072ce078c0'

  def test_hgrid_res0_5(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res0.5.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '2',
                         '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile), 'cda875e9d5f76eceeb98ac23d47574f1d8741be508a3f2d460b783b6ccd5c4f5'
  
  def test_hgrid_res0_5_equenh(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res0.5_equenh.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '2.0',
                         '--south_cutoff_row', '130', '--no_changing_meta',
                         '--write_subgrid_files', '--enhanced_equatorial'])
    assert sp.returncode == 0
    assert hashfile(outfile), 'e31104741e3d1180f1b14006b293892aeeedd94fda7afc57aeaa50dafebfe368'

  def test_hgrid_res0_25(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res0.25.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '4', '--rdp', '0.2',
                         '--south_cutoff_row', '83', '--write_subgrid_files', '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile), '5223dcca0e79f9d2a348a0c7b0946bab96f05e7669622ce7de1a5f570408a948'

  def test_hgrid_res0_125(self, tmpdir):
    outfile = tmpdir.join('ocean_hgrid_res0.125.nc')
    sp = subprocess.run([ogg_cmd,
                         '-f', outfile, '-r', '4', '--rdp', '0.2',
                         '--south_cutoff_row', '5', '--match_dy', '--even_j',
                         '--write_subgrid_files', '--no_changing_meta'])
    assert sp.returncode == 0
    assert hashfile(outfile), 'e3cb19861659332203c75ab7d12e9ce9150f89a7f33e578fbb7664f5a6b56156'

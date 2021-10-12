import pytest
import numpy
import ocean_grid_generator as ogg

def test_ogg_chksum(capsys):
    ogg.chksum(numpy.array([0, 1]), 'a')
    out, err = capsys.readouterr()
    assert 'fc62429c3e69001d65972cdeb94fb9aa18a7d9c16bc449e1e474e7e41bb95a7d          a min = 0.000000000000000 max = 1.000000000000000 mean = 0.500000000000000 sd = 0.500000000000000\n' == out
    pass




from .context import simulation
import numpy as np
import pytest


def test_ccd():
    """Tests that CCD is correctly generated"""
    with pytest.raises(TypeError):
        ccd = simulation.CCD()

    rows = 10
    columns = 20
    ccd = simulation.CCD(rows, columns)

    assert ccd.rows == rows
    assert ccd.columns == columns
    assert ccd.size() == (rows, columns)
    assert np.array_equal(ccd.image, np.zeros([rows, columns]))


def test_row_summation():
    rows = 2
    columns = 10
    ccd = simulation.CCD(rows, columns)
    ccd.set_image(np.ones([2, 10]))

    assert np.array_equal(ccd.shift_register, np.zeros([1, columns]))
    assert np.array_equal(ccd.shift_register[0, 0], 0)

    ccd.clock_row()
    assert np.array_equal(ccd.shift_register, np.ones([1, columns]))
    ccd.clock_row()
    assert np.array_equal(ccd.shift_register, np.ones([1, columns]) * 2)


def test_col_summation():
    rows = 1
    columns = 10
    ccd = simulation.CCD(rows, columns)
    ccd.set_image(np.ones([1, 10]))

    ccd.clock_row()

    ccd.clock_cloumn()
    assert ccd.well == 1
    test_array = np.ones([1, columns])
    test_array[0, -1] = 0
    assert np.array_equal(ccd.shift_register, test_array)

    ccd.clock_cloumn()
    assert ccd.well == 2
    test_array = np.ones([1, columns])
    test_array[0, -2:] = 0
    assert np.array_equal(ccd.shift_register, test_array)

    ccd.reset_shift_register()
    assert np.array_equal(ccd.shift_register, np.zeros([1, columns]))


def test_row_and_col_summation():
    rows = 2
    columns = 10
    ccd = simulation.CCD(rows, columns)
    ccd.set_image(np.ones([2, 10]))

    ccd.clock_row()
    ccd.clock_row()
    ccd.clock_cloumn()
    ccd.clock_cloumn()

    assert ccd.well == 4


def test_ccd_read():
    """tests reading of CCD"""
    rows = 30
    columns = 150
    ccd = simulation.CCD(rows, columns)
    with pytest.raises(ValueError):
        "test to ensure that image is same size as CCD"
        ccd.set_image(np.ones([1, 1]))

    ccd.set_image(np.ones([rows, columns]))

    assert np.array_equal(ccd.get_image(), np.ones([rows, columns]))


def test_read_well():
    """test that the well can be read"""
    ccd = simulation.CCD(512, 2048)
    ccd.well = np.ones([1, 1])
    assert ccd.read_well() == 1


def test_binning():
    """test binning functionality"""
    rows = 10
    columns = 20
    ccd = simulation.CCD(rows, columns)
    ccd.set_image(np.ones([rows, columns]))

    out_image = ccd.get_image(nrow=5, nrskip=0, nrbin=2, ncol=20, ncskip=0, ncbin=1)
    assert np.array_equal(out_image, np.ones([5, 20]) * 2)

    out_image = ccd.get_image(nrow=10, nrskip=0, nrbin=1, ncol=10, ncskip=0, ncbin=2)
    assert np.array_equal(out_image, np.ones([10, 10]) * 2)

    out_image = ccd.get_image(nrow=5, nrskip=0, nrbin=2, ncol=10, ncskip=0, ncbin=2)
    assert np.array_equal(out_image, np.ones([5, 10]) * 4)

    out_image = ccd.get_image(nrow=4, nrskip=0, nrbin=2, ncol=6, ncskip=0, ncbin=2)
    assert np.array_equal(out_image, np.ones([4, 6]) * 4)


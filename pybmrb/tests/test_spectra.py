#!/usr/bin/env python3

from pybmrb import Spectra
spect=Spectra()
def test_create_c13hsqc_peaklist():
    data = spect.create_c13hsqc_peaklist(bmrb_ids=27688)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_tocsy_peaklist():
    data = spect.create_tocsy_peaklist(bmrb_ids=30153)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_2d_peaklist():
    data = spect.create_2d_peaklist(bmrb_ids=25493, atom_x='H', atom_y='N')
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_n15hsqc_peaklist():
    data = spect.create_n15hsqc_peaklist(bmrb_ids=27688)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_n15hsqc():
    data = spect.n15hsqc(bmrb_ids=19152, show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_c13hsqc():
    data = spect.n15hsqc(bmrb_ids=18228, show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_tocsy():
    data = spect.tocsy(bmrb_ids=27435, show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_generic_2d():
    data = spect.generic_2d(bmrb_ids=18078, atom_x='H', atom_y='N', show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])
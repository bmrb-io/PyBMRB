#!/usr/bin/env python3

from pybmrb import Spectra


def test_create_c13hsqc_peaklist():
    data = Spectra.create_c13hsqc_peaklist(bmrb_ids=15060)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_tocsy_peaklist():
    data = Spectra.create_tocsy_peaklist(bmrb_ids=15060)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_2d_peaklist():
    data = Spectra.create_2d_peaklist(bmrb_ids=15060, atom_x='H', atom_y='N')
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_n15hsqc_peaklist():
    data = Spectra.create_n15hsqc_peaklist(bmrb_ids=15060)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_n15hsqc():
    data = Spectra.n15hsqc(bmrb_ids=15060,show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_c13hsqc():
    data = Spectra.n15hsqc(bmrb_ids=15060,show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_tocsy():
    data = Spectra.tocsy(bmrb_ids=15060,show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_generic_2d():
    data = Spectra.generic_2d(bmrb_ids=15060, atom_x='H', atom_y='N',show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])

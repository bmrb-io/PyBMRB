#!/usr/bin/env python3

from pybmrb import Spectra


def test_create_c13hsqc_peaklist():
    data = Spectra.create_c13hsqc_peaklist(bmrb_ids='27688')
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_tocsy_peaklist():
    data = Spectra.create_tocsy_peaklist(bmrb_ids='30153')
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_2d_peaklist():
    data = Spectra.create_2d_peaklist(bmrb_ids='25493', atom_x='H', atom_y='N')
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_create_n15hsqc_peaklist():
    data = Spectra.create_n15hsqc_peaklist(bmrb_ids='27688')
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_n15hsqc():
    data = Spectra.n15hsqc(bmrb_ids='19152', show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_c13hsqc():
    data = Spectra.n15hsqc(bmrb_ids='18228', show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_tocsy():
    data = Spectra.tocsy(bmrb_ids='27435', show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_generic_2d():
    data = Spectra.generic_2d(bmrb_ids='18078', atom_x='H', atom_y='N', show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 6
    for i in range(len(data) - 1):
        for j in range(len(data) - 1):
            if i != j:
                assert len(data[i]) == len(data[j])


def test_export_peak_list():
    data = Spectra.n15hsqc(bmrb_ids='17300', show_visualization=False)
    csv_dict = Spectra.export_peak_list(data, output_format='csv')
    assert type(csv_dict) is dict
    assert len(csv_dict.keys()) == 6
    n = len(csv_dict['sequence'])
    assert n > 0
    for k in csv_dict.keys():
        assert len(csv_dict[k]) == n
    sparky_dict = Spectra.export_peak_list(data, output_format='sparky')
    assert type(sparky_dict) is dict
    assert len(sparky_dict.keys()) == 3
    n = len(sparky_dict['Assignment'])
    assert n > 0
    for k in sparky_dict.keys():
        assert len(sparky_dict[k]) == n

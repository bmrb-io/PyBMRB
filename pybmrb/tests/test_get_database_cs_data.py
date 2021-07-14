#!/usr/bin/env python3

from pybmrb.get_database_cs_data import ChemicalShiftStatistics
import os
import ntpath

(script_path, script_name) = ntpath.split(os.path.realpath(__file__))
css = ChemicalShiftStatistics()

def test_get_data_from_api():
    dump = css._get_data_from_api('ALA', 'N')
    assert len(dump) == 2
    assert type(dump) is dict
    assert len(dump['columns']) == 13
    assert len(dump['data']) > 0


def test_get_data():
    dump = css.get_data(residue='ALA', atom='CB', filtered=True, ambiguity=1)
    assert type(dump) is tuple
    assert len(dump) == 2
    assert len(dump[0]) == 13
    assert len(dump[1]) > 0
    dump = css.get_data(residue='ALA', atom='*', filtered=True)
    assert type(dump) is tuple
    assert len(dump) == 2
    assert len(dump[0]) == 13
    assert len(dump[1]) > 0
    dump = css.get_data(residue='ALA', atom='H*', filtered=False)
    assert type(dump) is tuple
    assert len(dump) == 2
    assert len(dump[0]) == 13
    assert len(dump[1]) > 0


def test_get_data_from_bmrb():
    dump = css.get_data_from_bmrb(residue='ALA')
    assert type(dump) is tuple
    assert len(dump) == 2
    assert len(dump[0]) == 13
    assert len(dump[1]) > 0


def test_list_do_dict():
    dump = css.get_data_from_bmrb(residue='ALA')
    dump2 = css._list_to_dict(dump[0], dump[1])
    assert type(dump2) is dict


def test_get_2d_chemical_shifts():
    dump = css.get_2d_chemical_shifts(residue='CYS', atom1='N', atom2='CB')
    assert type(dump) is tuple
    assert len(dump[0]) == len(dump[1])


def test_get_filtered_data_from_bmrb():
    dump = css.get_filtered_data_from_bmrb(residue='CYS', atom='CB',
                                                               filtering_rules=[('CA', 64.5), ('H', 7.6)])
    assert type(dump) is list
    assert len(dump) > 0


def test_get_statistics():
    dump = css.get_statistics(residue='TYR', atom='CB')
    assert type(dump) is dict
    for key in dump:
        assert len(dump[key].keys()) == 7

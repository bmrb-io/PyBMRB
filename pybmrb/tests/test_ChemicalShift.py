#!/usr/bin/env python3

import pynmrstar
from pybmrb import ChemicalShift
import os
import ntpath

(script_path, script_name) = ntpath.split(os.path.realpath(__file__))


def test__from_pynmrstar_entry_object():
    entry_obj1 = pynmrstar.Entry.from_file('{}/test_data/MyData.str'.format(script_path))
    entry_obj2 = pynmrstar.Entry.from_database('19798')
    entry_obj3 = pynmrstar.Entry.from_database('18857')
    assert len(ChemicalShift._from_pynmrstar_entry_object(entry_obj1, data_set_id='MyData')) == 1
    assert type(ChemicalShift._from_pynmrstar_entry_object(entry_obj1, data_set_id='MyData')) is dict
    assert len(ChemicalShift._from_pynmrstar_entry_object(entry_obj2, data_set_id='19798')) == 1
    assert type(ChemicalShift._from_pynmrstar_entry_object(entry_obj2, data_set_id='19798')) is dict
    assert len(ChemicalShift._from_pynmrstar_entry_object(entry_obj3, data_set_id='18857')) == 33
    assert type(ChemicalShift._from_pynmrstar_entry_object(entry_obj3, data_set_id='18857')) is dict
    assert len(
        ChemicalShift._from_pynmrstar_entry_object(entry_obj1, data_set_id='MyData', auth_tag=True)) == 1
    assert type(
        ChemicalShift._from_pynmrstar_entry_object(entry_obj1, data_set_id='MyData', auth_tag=True)) is dict
    assert len(ChemicalShift._from_pynmrstar_entry_object(entry_obj2, data_set_id='19798', auth_tag=True)) == 1
    assert type(
        ChemicalShift._from_pynmrstar_entry_object(entry_obj2, data_set_id='19798', auth_tag=True)) is dict
    assert len(
        ChemicalShift._from_pynmrstar_entry_object(entry_obj3, data_set_id='18857', auth_tag=True)) == 33
    assert type(
        ChemicalShift._from_pynmrstar_entry_object(entry_obj3, data_set_id='18857', auth_tag=True)) is dict


def test_from_entry_object():
    ent1 = pynmrstar.Entry.from_database('15000')
    ent2 = pynmrstar.Entry.from_database('18857')
    cs_data1 = ChemicalShift.from_entry_object(entry_objects=ent1, data_set_id=None)
    cs_data2 = ChemicalShift.from_entry_object(entry_objects=ent2)
    assert len(cs_data1) == 1
    assert type(cs_data1) is dict
    assert len(cs_data2) == 33
    assert type(cs_data2) is dict


def test_from_file():
    cs_data1 = ChemicalShift.from_file(input_file_names='{}/test_data/MyData.str'.format(script_path))
    cs_data2 = ChemicalShift.from_file(input_file_names='{}/test_data/bmr18857_3.str'.format(script_path))
    assert len(cs_data1) == 1
    assert type(cs_data1) is dict
    assert len(cs_data2) == 33
    assert type(cs_data2) is dict
    cs_data3 = ChemicalShift.from_file(input_file_names=['{}/test_data/MyData.str'.format(script_path),
                                                         '{}/test_data/bmr18857_3.str'.format(script_path)])
    assert len(cs_data3) == 34
    assert type(cs_data3) is dict
    cs_data1 = ChemicalShift.from_file(input_file_names='{}/test_data/MyData.str'.format(script_path), auth_tag=True)
    cs_data2 = ChemicalShift.from_file(input_file_names='{}/test_data/bmr18857_3.str'.format(script_path),
                                       auth_tag=True)
    assert len(cs_data1) == 1
    assert type(cs_data1) is dict
    assert len(cs_data2) == 33
    assert type(cs_data2) is dict
    cs_data3 = ChemicalShift.from_file(input_file_names=['{}/test_data/MyData.str'.format(script_path),
                                                         '{}/test_data/bmr18857_3.str'.format(script_path)],
                                       auth_tag=True)
    assert len(cs_data3) == 34
    assert type(cs_data3) is dict


def test_from_bmrb():
    cs_data1 = ChemicalShift.from_bmrb(bmrb_ids='15060', auth_tag=True)
    cs_data2 = ChemicalShift.from_bmrb(bmrb_ids='18857', auth_tag=True)
    assert len(cs_data1) == 1
    assert type(cs_data1) is dict
    assert len(cs_data2) == 33
    assert type(cs_data2) is dict
    cs_data3 = ChemicalShift.from_bmrb(bmrb_ids=['15000', 18857], auth_tag=True)
    assert len(cs_data3) == 34
    assert type(cs_data3) is dict
    cs_data1 = ChemicalShift.from_bmrb(bmrb_ids='50957')
    cs_data2 = ChemicalShift.from_bmrb(bmrb_ids='18857')
    assert len(cs_data1) == 1
    assert type(cs_data1) is dict
    assert len(cs_data2) == 33
    assert type(cs_data2) is dict
    cs_data3 = ChemicalShift.from_bmrb(bmrb_ids=['28038', 18857])
    assert len(cs_data3) == 34
    assert type(cs_data3) is dict

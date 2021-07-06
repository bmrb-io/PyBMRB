import pytest
import pynmrstar
from pybmrb2.get_entry_cs_data import ChemicalShift
import os
import ntpath

(script_path, script_name) = ntpath.split(os.path.realpath(__file__))



@pytest.fixture
def get_entry_object():
        return [pynmrstar.Entry.from_file('{}/test_data/MyData.str'.format(script_path).format(script_path)),
                pynmrstar.Entry.from_database('15060'),
                pynmrstar.Entry.from_database('18857')]




def test_from_entry(get_entry_object):
    assert len(ChemicalShift.from_entry(get_entry_object[0],data_set_id='MyData')) == 1
    assert type(ChemicalShift.from_entry(get_entry_object[0],data_set_id='MyData')) is dict
    assert len(ChemicalShift.from_entry(get_entry_object[1], data_set_id='15560')) == 1
    assert type(ChemicalShift.from_entry(get_entry_object[1], data_set_id='15560')) is dict
    assert len(ChemicalShift.from_entry(get_entry_object[2], data_set_id='18857')) == 33
    assert type(ChemicalShift.from_entry(get_entry_object[2], data_set_id='18857')) is dict
    assert len(ChemicalShift.from_entry(get_entry_object[0], data_set_id='MyData',auth_tag=True)) == 1
    assert type(ChemicalShift.from_entry(get_entry_object[0], data_set_id='MyData',auth_tag=True)) is dict
    assert len(ChemicalShift.from_entry(get_entry_object[1], data_set_id='15560',auth_tag=True)) == 1
    assert type(ChemicalShift.from_entry(get_entry_object[1], data_set_id='15560',auth_tag=True)) is dict
    assert len(ChemicalShift.from_entry(get_entry_object[2], data_set_id='18857',auth_tag=True)) == 33
    assert type(ChemicalShift.from_entry(get_entry_object[2], data_set_id='18857',auth_tag=True)) is dict



def test_from_file():
    assert len(ChemicalShift.from_file('{}/test_data/MyData.str'.format(script_path).format(script_path))) == 1
    assert type(ChemicalShift.from_file('{}/test_data/MyData.str'.format(script_path))) is dict
    assert len(ChemicalShift.from_file('{}/test_data/bmr18857_3.str'.format(script_path))) == 33
    assert type(ChemicalShift.from_file('{}/test_data/bmr18857_3.str'.format(script_path))) is dict
    assert len(ChemicalShift.from_file(['{}/test_data/MyData.str'.format(script_path), '{}/test_data/bmr18857_3.str'.format(script_path)])) == 34
    assert type(ChemicalShift.from_file(['{}/test_data/MyData.str'.format(script_path), '{}/test_data/bmr18857_3.str'.format(script_path)])) is dict
    assert len(ChemicalShift.from_file('{}/test_data/MyData.str'.format(script_path), auth_tag=True)) == 1
    assert type(ChemicalShift.from_file('{}/test_data/MyData.str'.format(script_path), auth_tag=True)) is dict
    assert len(ChemicalShift.from_file('{}/test_data/bmr18857_3.str'.format(script_path), auth_tag=True)) == 33
    assert type(ChemicalShift.from_file('{}/test_data/bmr18857_3.str'.format(script_path), auth_tag=True)) is dict
    assert len(ChemicalShift.from_file(['{}/test_data/MyData.str'.format(script_path), '{}/test_data/bmr18857_3.str'.format(script_path)], auth_tag=True)) == 34
    assert type(ChemicalShift.from_file(['{}/test_data/MyData.str'.format(script_path), '{}/test_data/bmr18857_3.str'.format(script_path)], auth_tag=True)) is dict


def test_from_bmrb():
    assert len(ChemicalShift.from_bmrb('15060',auth_tag=True)) == 1
    assert type(ChemicalShift.from_bmrb(15060,auth_tag=True)) is dict
    assert len(ChemicalShift.from_bmrb('18857',auth_tag=True)) == 33
    assert type(ChemicalShift.from_bmrb(18857,auth_tag=True)) is dict
    assert len(ChemicalShift.from_bmrb(['15060',18857],auth_tag=True)) == 34
    assert type(ChemicalShift.from_bmrb([15060,18857],auth_tag=True)) is dict
    assert len(ChemicalShift.from_bmrb('15060')) == 1
    assert type(ChemicalShift.from_bmrb(15060)) is dict
    assert len(ChemicalShift.from_bmrb('18857')) == 33
    assert type(ChemicalShift.from_bmrb(18857)) is dict
    assert len(ChemicalShift.from_bmrb(['15060', 18857])) == 34
    assert type(ChemicalShift.from_bmrb([15060, 18857])) is dict

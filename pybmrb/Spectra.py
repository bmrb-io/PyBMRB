#!/usr/bin/env python3
"""
This module is used to visualize the one dimensional chemical shift list from BMRB entry or NMR-STAR file as a 
two dimensional NMR spectrum. It simulates peak positions of \u00b9\u2075N-HSQC, \u00b9\u00b3C-HSQC and
\u00b9H-\u00b9H-TOCSY. It can also simulate a generic 2D spectrum between any two given pair of atoms. This module is 
useful to compare user data with any BMRB entry/list of entries as a overlaid NMR spectra.
"""

import logging
import csv

import pynmrstar

from pybmrb import ChemicalShift, ChemicalShiftStatistics
import plotly.express as px
from typing import Union, List, Optional


def create_c13hsqc_peaklist(bmrb_ids: Optional[Union[str, List[str], int, List[int]]] = None,
                            entry_objects: Optional[Union[pynmrstar.Entry, List[pynmrstar.Entry]]] = None,
                            input_file_names: Optional[Union[str, List[str]]] = None,
                            auth_tag: Optional[bool] = False,
                            draw_trace: Optional[bool] = False) -> tuple:
    """
    Converts one dimensional chemical shifts from BMRB entries/NMR-STAR files/PyNMRSTAR entry objects into
    \u00b9\u00b3C-HSQC peak list

    :param bmrb_ids: single BMRB entry ID or multiple BMRB entry IDs as list, defaults to None
    :type bmrb_ids: str/int/list, optional
    :param input_file_names: single input file name  or multiple input file names as list, defaults to None
    :type input_file_names: str, optional
    :param entry_objects: single PyNMRSTAR entry object  or multiple PyNMRSTAR entry objects as list, defaults to None
    :type entry_objects: PyNMRSTAR object/list, optional
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :param draw_trace: draw line connecting peaks from residues at the same sequence location in different
        data sets, defaults to False
    :type draw_trace: bool optional
    :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track)
        if draw_trace is True cs_track={ matching atoms:[cs_values]} else cs_track={}
    :rtype: tuple
    """
    if bmrb_ids is None and input_file_names is None and entry_objects is None:
        logging.error('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
        raise TypeError('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
    ch_atoms = {'ALA': [('HA', 'CA'), ('HB*', 'CB')],
                'ARG': [('HA', 'CA'), ('HB*', 'CB'), ('HG*', 'CG'), ('HD*', 'CD')],
                'ASN': [('HA', 'CA'), ('HB*', 'CB')],
                'ASP': [('HA', 'CA'), ('HB*', 'CB')],
                'CYS': [('HA', 'CA'), ('HB*', 'CB')],
                'GLN': [('HA', 'CA'), ('HB*', 'CB'), ('HG*', 'CG')],
                'GLU': [('HA', 'CA'), ('HB*', 'CB'), ('HG*', 'CG')],
                'GLY': [('HA*', 'CA')],
                'HIS': [('HA', 'CA'), ('HB*', 'CB'), ('HD2', 'CD2'), ('HE1', 'CE1')],
                'ILE': [('HA', 'CA'), ('HB', 'CB'), ('HG2*', 'CG2'), ('HG1*', 'CG1'), ('HD1*', 'CD1')],
                'LEU': [('HA', 'CA'), ('HB*', 'CB'), ('HG', 'CG'), ('HD1*', 'CD1'), ('HD2*', 'CD2')],
                'LYS': [('HA', 'CA'), ('HB*', 'CB'), ('HG*', 'CG'), ('HD*', 'CD'), ('HE*', 'CE')],
                'MET': [('HA', 'CA'), ('HB*', 'CB'), ('HG*', 'CG'), ('HE*', 'CE')],
                'PHE': [('HA', 'CA'), ('HB*', 'CB'), ('HD1', 'CD1'), ('HE1', 'CE1'), ('HD2', 'CD2'), ('HE2', 'CE2'),
                        ('HZ', 'CZ')],
                'PRO': [('HA', 'CA'), ('HB*', 'CB'), ('HG*', 'CG'), ('HD*', 'CD')],
                'SER': [('HA', 'CA'), ('HB*', 'CB')],
                'TRP': [('HA', 'CA'), ('HB*', 'CB'), ('HD1', 'CD1'), ('HE3', 'CE3'), ('HZ2', 'CZ2'),
                        ('HH2', 'CH2')],
                'TYR': [('HA', 'CA'), ('HB*', 'CB'), ('HD1', 'CD1'), ('HE1', 'CE1'), ('HD2', 'CD2'),
                        ('HE2', 'CE2')],
                'THR': [('HA', 'CA'), ('HB', 'CB'), ('HG2*', 'CG2')],
                'VAL': [('HA', 'CA'), ('HB', 'CB'), ('HG1*', 'CG1'), ('HG2*', 'CG2')],
                }
    cs_data = {}
    cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids=bmrb_ids, auth_tag=auth_tag)
    cs_data.update(cs_data_bmrb)
    if input_file_names is not None:
        cs_data_file = ChemicalShift.from_file(input_file_names=input_file_names, auth_tag=auth_tag)
        cs_data.update(cs_data_file)
    if entry_objects is not None:
        cs_data_obj = ChemicalShift.from_entry_object(entry_objects=entry_objects, auth_tag=auth_tag)
        cs_data.update(cs_data_obj)
    data_set = []
    x = []
    y = []
    res = []
    info = []
    atom_ids = {}
    for data_id in cs_data.keys():
        for chain in cs_data[data_id].keys():
            for seq_no in cs_data[data_id][chain]['seq_ids']:
                residue = cs_data[data_id][chain][seq_no][list(cs_data[data_id][chain][seq_no].keys())[0]][0]
                ch_list = ch_atoms[residue]
                for ch_atom in ch_list:
                    for atm_x in cs_data[data_id][chain][seq_no].keys():
                        if ('*' in ch_atom[0]
                            and ch_atom[0][:-1] == atm_x[:len(ch_atom[0][:-1])]) \
                                or (ch_atom[0] == atm_x):
                            cs_x = cs_data[data_id][chain][seq_no][atm_x][2]
                            for atm_y in cs_data[data_id][chain][seq_no].keys():
                                if ('*' in ch_atom[1]
                                    and ch_atom[1][:-1] == atm_y[:len(ch_atom[1][:-1])]) \
                                        or (ch_atom[1] == atm_y):
                                    cs_y = cs_data[data_id][chain][seq_no][atm_y][2]
                                    data_set.append(data_id)
                                    x.append(cs_x)
                                    y.append(cs_y)
                                    res.append(residue)
                                    tag = '{}-{}-{}-{}-{}-{}'.format(data_id, chain, seq_no, residue, atm_x, atm_y)
                                    info.append(tag)
                                    atom_id = '{}-{}-{}-{}-{}'.format(chain, seq_no, residue, atm_x, atm_y)
                                    if draw_trace:
                                        if atom_id not in atom_ids.keys():
                                            atom_ids[atom_id] = [[], []]
                                        atom_ids[atom_id][0].append(cs_x)
                                        atom_ids[atom_id][1].append(cs_y)
    cs_track = {}
    if draw_trace:
        for k in atom_ids.keys():
            if len(atom_ids[k][0]) > 1:
                cs_track[k] = atom_ids[k]
    return x, y, data_set, info, res, cs_track


def create_tocsy_peaklist(bmrb_ids: Optional[Union[str, List[str], int, List[int]]] = None,
                          entry_objects: Optional[Union[pynmrstar.Entry, List[pynmrstar.Entry]]] = None,
                          input_file_names: Optional[Union[str, List[str]]] = None,
                          auth_tag: Optional[bool] = False,
                          draw_trace: Optional[bool] = False) -> tuple:
    """
    Converts one dimensional chemical shifts from BMRB entries/NMR-STAR files/PyNMRSTAR entry objects
    into  into \u00b9H-\u00b9H-TOCSY peak list

    :param bmrb_ids: single BMRB entry ID or multiple BMRB entry IDs as list, defaults to None
    :type bmrb_ids: str/int/list, optional
    :param input_file_names: single input file name  or multiple input file names as list, defaults to None
    :type input_file_names: str, optional
    :param entry_objects: single PyNMRSTAR entry object  or multiple PyNMRSTAR entry objects as list, defaults to None
    :type entry_objects: PyNMRSTAR object/list, optional
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :param draw_trace: draw line connecting peaks from residues at the same sequence location in different
        data sets, defaults to False
    :type draw_trace: bool optional
    :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track)
        if draw_trace is True cs_track={ matching atoms:[cs_values]} else cs_track={}
    :rtype: tuple
    """
    if bmrb_ids is None and input_file_names is None and entry_objects is None:
        logging.error('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
        raise TypeError('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
    cs_data = {}
    cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids=bmrb_ids, auth_tag=auth_tag)
    cs_data.update(cs_data_bmrb)
    if input_file_names is not None:
        cs_data_file = ChemicalShift.from_file(input_file_names=input_file_names, auth_tag=auth_tag)
        cs_data.update(cs_data_file)
    if entry_objects is not None:
        cs_data_obj = ChemicalShift.from_entry_object(entry_objects=entry_objects, auth_tag=auth_tag)
        cs_data.update(cs_data_obj)
    data_set = []
    x = []
    y = []
    res = []
    info = []
    atom_ids = {}
    for data_id in cs_data.keys():
        for chain in cs_data[data_id].keys():
            for seq_no in cs_data[data_id][chain]['seq_ids']:
                residue = cs_data[data_id][chain][seq_no][list(cs_data[data_id][chain][seq_no].keys())[0]][0]
                ch_list = [('H*', 'H*')]
                for ch_atom in ch_list:
                    for atm_x in cs_data[data_id][chain][seq_no].keys():
                        if ('*' in ch_atom[0]
                            and ch_atom[0][:-1] == atm_x[:len(ch_atom[0][:-1])]) \
                                or (ch_atom[0] == atm_x):
                            cs_x = cs_data[data_id][chain][seq_no][atm_x][2]
                            for atm_y in cs_data[data_id][chain][seq_no].keys():
                                if ('*' in ch_atom[1]
                                    and ch_atom[1][:-1] == atm_y[:len(ch_atom[1][:-1])]) \
                                        or (ch_atom[1] == atm_y):
                                    cs_y = cs_data[data_id][chain][seq_no][atm_y][2]
                                    data_set.append(data_id)
                                    x.append(cs_x)
                                    y.append(cs_y)
                                    res.append(residue)
                                    tag = '{}-{}-{}-{}-{}-{}'.format(data_id, chain, seq_no, residue, atm_x, atm_y)
                                    info.append(tag)
                                    atom_id = '{}-{}-{}-{}-{}'.format(chain, seq_no, residue, atm_x, atm_y)
                                    if draw_trace:
                                        if atom_id not in atom_ids.keys():
                                            atom_ids[atom_id] = [[], []]
                                        atom_ids[atom_id][0].append(cs_x)
                                        atom_ids[atom_id][1].append(cs_y)
    cs_track = {}
    if draw_trace:
        for k in atom_ids.keys():
            if len(atom_ids[k][0]) > 1:
                cs_track[k] = atom_ids[k]
    return x, y, data_set, info, res, cs_track


def create_2d_peaklist(atom_x: str,
                       atom_y: str,
                       bmrb_ids: Optional[Union[str, List[str], int, List[int]]] = None,
                       input_file_names: Optional[Union[str, List[str]]] = None,
                       entry_objects: Optional[Union[pynmrstar.Entry, List[pynmrstar.Entry]]] = None,
                       auth_tag: Optional[bool] = False,
                       draw_trace: Optional[bool] = False,
                       include_preceding: Optional[bool] = False,
                       include_next: Optional[bool] = False,
                       legend: Optional[str] = None) -> tuple:
    """
     Converts one dimensional chemical shifts from BMRB entries/NMR-STAR files/PyNMRSTAR entry objects
     into  into generic 2D peak list

    :param atom_x: atom name for X coordinate in IUPAC format
    :type atom_x: str
    :param atom_y: atom name for Y coordinate in IUPAC format
    :type atom_y: str
    :param bmrb_ids: single BMRB entry ID or multiple BMRB entry IDs as list, defaults to None
    :type bmrb_ids: str/int/list, optional
    :param input_file_names: single input file name  or multiple input file names as list, defaults to None
    :type input_file_names: str, optional
    :param entry_objects: single PyNMRSTAR entry object  or multiple PyNMRSTAR entry objects as list, defaults to None
    :type entry_objects: PyNMRSTAR object/list, optional
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :param draw_trace: draw line connecting peaks from residues at the same sequence location in different
        data sets, defaults to False
    :type draw_trace: bool optional
    :param include_preceding: include peaks from i-1 residue on the Y axis, defaults to False
    :type include_preceding: bool, optional
    :param include_next: include peaks from i+i residue on the Y axis, defaults to False
    :type include_next: bool, optional
    :param legend: legends are disabled by default. Residue types are indicated by  color and  data sets are
        indicated by symbol, displaying the combination of both will create a very long list of legend. Optionally
        either 'residue' or 'dataset' can be used to color code the scatter plot by residue type
        or data set and display the legend, defaults to None
    :type legend: str, optional
    :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track,psn,seq_track)
        if draw_trace is True cs_track={ matching atoms:[cs_values]} else cs_track={}
    :rtype: tuple
    """
    if bmrb_ids is None and input_file_names is None and entry_objects is None:
        logging.error('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
        raise TypeError('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')

    cs_data = {}
    cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids, auth_tag=auth_tag)
    cs_data.update(cs_data_bmrb)
    if input_file_names is not None:
        cs_data_file = ChemicalShift.from_file(input_file_names, auth_tag=auth_tag)
        cs_data.update(cs_data_file)
    if entry_objects is not None:
        cs_data_obj = ChemicalShift.from_entry_object(entry_objects=entry_objects, auth_tag=auth_tag)
        cs_data.update(cs_data_obj)
    data_set = []
    x = []
    y = []
    res = []
    info = []
    psn = []
    atom_ids = {}
    seq_trace = [[], [], []]
    for data_id in cs_data.keys():
        for chain in cs_data[data_id].keys():
            for seq_no in cs_data[data_id][chain]['seq_ids']:
                residue = cs_data[data_id][chain][seq_no][list(cs_data[data_id][chain][seq_no].keys())[0]][0]
                ch_list = [(atom_x, atom_y)]
                for ch_atom in ch_list:
                    for atm_x in cs_data[data_id][chain][seq_no].keys():
                        if ('*' in ch_atom[0]
                            and ch_atom[0][:-1] == atm_x[:len(ch_atom[0][:-1])]) \
                                or (ch_atom[0] == atm_x):
                            cs_x = cs_data[data_id][chain][seq_no][atm_x][2]
                            for atm_y in cs_data[data_id][chain][seq_no].keys():
                                if ('*' in ch_atom[1]
                                    and ch_atom[1][:-1] == atm_y[:len(ch_atom[1][:-1])]) \
                                        or (ch_atom[1] == atm_y):
                                    if include_preceding:
                                        try:
                                            residue_p = cs_data[data_id][chain][seq_no - 1][
                                                list(cs_data[data_id][chain][seq_no - 1].keys())[0]][0]
                                            cs_y = cs_data[data_id][chain][seq_no - 1][atm_y][2]
                                            data_set.append(data_id)
                                            x.append(cs_x)
                                            y.append(cs_y)
                                            seq_trace[0].append(cs_x)
                                            seq_trace[1].append(cs_y)
                                            seq_trace[2].append(seq_no - 1)
                                            res.append(residue)
                                            if legend == 'psn':
                                                psn.append('Preceding')
                                            else:
                                                psn.append(0.5)
                                            tag = '{}-{}-{}-{}-{}-{}-{}-{}'.format(data_id, chain, seq_no, residue,
                                                                                   seq_no - 1, residue_p, atm_x,
                                                                                   atm_y)
                                            info.append(tag)
                                            atom_id = '{}-{}-{}-{}-{}-{}-{}'.format(chain, seq_no, residue, seq_no - 1,
                                                                                    residue_p, atm_x, atm_y)
                                            if draw_trace:
                                                if atom_id not in atom_ids.keys():
                                                    atom_ids[atom_id] = [[], []]
                                                atom_ids[atom_id][0].append(cs_x)
                                                atom_ids[atom_id][1].append(cs_y)
                                        except KeyError:
                                            pass
                                    cs_y = cs_data[data_id][chain][seq_no][atm_y][2]
                                    data_set.append(data_id)
                                    x.append(cs_x)
                                    y.append(cs_y)
                                    seq_trace[0].append(cs_x)
                                    seq_trace[1].append(cs_y)
                                    seq_trace[2].append(seq_no)
                                    res.append(residue)
                                    if legend == 'psn':
                                        psn.append('Same')
                                    else:
                                        psn.append(1)
                                    tag = '{}-{}-{}-{}-{}-{}'.format(data_id, chain, seq_no, residue, atm_x, atm_y)
                                    info.append(tag)
                                    atom_id = '{}-{}-{}-{}-{}'.format(chain, seq_no, residue, atm_x, atm_y)
                                    if draw_trace:
                                        if atom_id not in atom_ids.keys():
                                            atom_ids[atom_id] = [[], []]
                                        atom_ids[atom_id][0].append(cs_x)
                                        atom_ids[atom_id][1].append(cs_y)
                                    if include_next:
                                        try:
                                            residue_p = cs_data[data_id][chain][seq_no + 1][
                                                list(cs_data[data_id][chain][seq_no + 1].keys())[0]][0]
                                            cs_y = cs_data[data_id][chain][seq_no + 1][atm_y][2]
                                            data_set.append(data_id)
                                            x.append(cs_x)
                                            y.append(cs_y)
                                            seq_trace[0].append(cs_x)
                                            seq_trace[1].append(cs_y)
                                            seq_trace[2].append(seq_no + 1)
                                            res.append(residue)
                                            if legend == 'psn':
                                                psn.append('Next')
                                            else:
                                                psn.append(0.5)
                                            tag = '{}-{}-{}-{}-{}-{}-{}-{}'.format(data_id, chain, seq_no, residue,
                                                                                   seq_no + 1, residue_p, atm_x,
                                                                                   atm_y)
                                            info.append(tag)
                                            atom_id = '{}-{}-{}-{}-{}-{}-{}'.format(chain, seq_no, residue, seq_no + 1,
                                                                                    residue_p, atm_x, atm_y)
                                            if draw_trace:
                                                if atom_id not in atom_ids.keys():
                                                    atom_ids[atom_id] = [[], []]
                                                atom_ids[atom_id][0].append(cs_x)
                                                atom_ids[atom_id][1].append(cs_y)
                                        except KeyError:
                                            pass

    cs_track = {}
    if draw_trace:
        for k in atom_ids.keys():
            if len(atom_ids[k][0]) > 1:
                cs_track[k] = atom_ids[k]
    return x, y, data_set, info, res, cs_track, psn, seq_trace


def create_n15hsqc_peaklist(bmrb_ids: Optional[Union[str, List[str], int, List[int]]] = None,
                            input_file_names: Optional[Union[str, List[str]]] = None,
                            entry_objects: Optional[Union[pynmrstar.Entry, List[pynmrstar.Entry]]] = None,
                            auth_tag: Optional[bool] = False,
                            draw_trace: Optional[bool] = False,
                            include_sidechain: Optional[bool] = True) -> tuple:
    """
    Converts one dimensional chemical shifts from BMRB entries/NMR-STAR files/PyNMRSTAR
    entry objects into  \u00b9\u2075N-HSQC peak list

    :param bmrb_ids: single BMRB entry ID or multiple BMRB entry IDs as list, defaults to None
    :type bmrb_ids: str/int/list, optional
    :param input_file_names: single input file name  or multiple input file names as list, defaults to None
    :type input_file_names: str, optional
    :param entry_objects: single PyNMRSTAR entry object  or multiple PyNMRSTAR entry objects as list, defaults to None
    :type entry_objects: PyNMRSTAR object/list, optional
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :param draw_trace: draw line connecting peaks from residues at the same sequence location in different
        data sets, defaults to False
    :type draw_trace: bool optional
    :param include_sidechain: include peaks from side chains, defaults to True
    :type include_sidechain: bool, optional
    :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track)
        if draw_trace is True cs_track={ matching atoms:[cs_values]} else cs_track={}
    :rtype: tuple
    """
    if bmrb_ids is None and input_file_names is None and entry_objects is None:
        logging.error('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
        raise TypeError('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
    atom_x = 'H'
    atom_y = 'N'
    sidechain_nh_atoms = {'ARG': {
        'ARG-HH11-NH1': ['HH11', 'NH1'],
        'ARG-HH12-NH1': ['HH12', 'NH1'],
        'ARG-HH21-NH2': ['HH21', 'NH2'],
        'ARG-HH22-NH2': ['HH22', 'NH2'],
        'ARG-HE-NE': ['HE', 'NE']},
        'GLN': {
            'GLN-HE21-NE2': ['HE21', 'NE2'],
            'GLN-HE22-NE2': ['HE22', 'NE2']},
        'ASN': {
            'ASN-HD21-ND2': ['HD21', 'ND2'],
            'ASN-HD22-ND2': ['HD22', 'ND2']},
        'HIS': {
            'HIS-HD1-ND1': ['HD1', 'ND1'],
            'HIS-HE2-ND1': ['HE2', 'NE2']},
        'TRP': {
            'TRP-HE1-NE1': ['HE1', 'NE1']},
        'LYS': {
            'LYS-HZ-NZ': ['HZ', 'NZ'],
            'LYS-HZ1-NZ': ['HZ1', 'NZ'],
            'LYS-HZ2-NZ': ['HZ2', 'NZ'],
            'LYS-HZ3-NZ': ['HZ3', 'NZ']}
    }
    cs_data = {}
    if bmrb_ids is not None:
        cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids=bmrb_ids, auth_tag=auth_tag)
        cs_data.update(cs_data_bmrb)
    if input_file_names is not None:
        cs_data_file = ChemicalShift.from_file(input_file_names=input_file_names, auth_tag=auth_tag)
        cs_data.update(cs_data_file)
    if entry_objects is not None:
        cs_data_obj = ChemicalShift.from_entry_object(entry_objects=entry_objects, auth_tag=auth_tag)
        cs_data.update(cs_data_obj)
    data_set = []
    x = []
    y = []
    res = []
    info = []
    atom_ids = {}
    for data_id in cs_data.keys():
        for chain in cs_data[data_id].keys():
            for seq_no in cs_data[data_id][chain]['seq_ids']:
                try:
                    x_cs = cs_data[data_id][chain][seq_no][atom_x][2]
                    y_cs = cs_data[data_id][chain][seq_no][atom_y][2]
                    residue = cs_data[data_id][chain][seq_no][atom_y][0]
                    res.append(residue)
                    tag = '{}-{}-{}-{}-H-N'.format(data_id, chain, seq_no, residue)
                    data_set.append(data_id)
                    atom_id = '{}-{}-{}'.format(chain, seq_no, residue)
                    if draw_trace:
                        if atom_id not in atom_ids.keys():
                            atom_ids[atom_id] = [[], []]
                        atom_ids[atom_id][0].append(x_cs)
                        atom_ids[atom_id][1].append(y_cs)
                    x.append(x_cs)
                    y.append(y_cs)
                    info.append(tag)
                    if include_sidechain:
                        if residue in sidechain_nh_atoms.keys():
                            for atom_list in sidechain_nh_atoms[residue]:
                                ax = sidechain_nh_atoms[residue][atom_list][0]
                                ay = sidechain_nh_atoms[residue][atom_list][1]
                                try:
                                    x_cs = cs_data[data_id][chain][seq_no][ax][2]
                                    y_cs = cs_data[data_id][chain][seq_no][ay][2]
                                    res.append(residue)
                                    tag = '{}-{}-{}-{}'.format(data_id, chain, seq_no,
                                                               atom_list)
                                    data_set.append(data_id)
                                    atom_id = '{}-{}-{}'.format(chain, seq_no,
                                                                atom_list)
                                    if draw_trace:
                                        if atom_id not in atom_ids.keys():
                                            atom_ids[atom_id] = [[], []]
                                        atom_ids[atom_id][0].append(x_cs)
                                        atom_ids[atom_id][1].append(y_cs)
                                    x.append(x_cs)
                                    y.append(y_cs)
                                    info.append(tag)
                                except KeyError:
                                    logging.debug('Data not found:{},{},{}'.format(data_id, chain, seq_no))
                except KeyError:
                    logging.debug('Data not found:{},{},{}'.format(data_id, chain, seq_no))
    cs_track = {}
    if draw_trace:
        for k in atom_ids.keys():
            if len(atom_ids[k][0]) > 1:
                cs_track[k] = atom_ids[k]
    return x, y, data_set, info, res, cs_track


def n15hsqc(bmrb_ids: Optional[Union[str, List[str], int, List[int]]] = None,
            input_file_names: Optional[Union[str, List[str]]] = None,
            entry_objects: Optional[Union[pynmrstar.Entry, List[pynmrstar.Entry]]] = None,
            auth_tag: Optional[bool] = False,
            legend: Optional[str] = None,
            draw_trace: Optional[bool] = False,
            include_sidechain: Optional[bool] = True,
            peak_list: Optional[str] = None,
            output_format: Optional[str] = 'html',
            output_file: Optional[str] = None,
            output_image_width: Optional[int] = 800,
            output_image_height: Optional[int] = 600,
            show_visualization: Optional[bool] = True) -> tuple:
    """
    Plots \u00b9\u2075N-HSQC spectrum  for a given  BMRB entry/NMR-STAR file/PyNMRSTAR entry object;
    This function can be used to compare different data sets as overlaid NMR Spectra. It overlays \u00b9\u2075N-HSQC
    for a given list of BMRB entries/NMR-STAR files/PyNMRSTAR entry objects and draw lines connecting peaks
    from residues at the same sequence location in different data sets

    :param bmrb_ids: single BMRB entry ID or multiple BMRB entry IDs as list, defaults to None
    :type bmrb_ids: str/int/list, optional
    :param input_file_names: single input file name  or multiple input file names as list, defaults to None
    :type input_file_names: str, optional
    :param entry_objects: single PyNMRSTAR entry object  or multiple PyNMRSTAR entry objects as list, defaults to None
    :type entry_objects: PyNMRSTAR object/list, optional
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :param legend: legends are disabled by default. Residue types are indicated by  color and  data sets are
        indicated by symbol, displaying the combination of both will create a very long list of legend. Optionally
        either 'residue' or 'dataset' can be used to color code the scatter plot by residue type
        or data set and display the legend, defaults to None
    :type legend: str, optional
    :param draw_trace: draw line connecting peaks from residues at the same sequence location in different
        data sets, defaults to False
    :type draw_trace: bool optional
    :param include_sidechain: include peaks from side chains, defaults to True
    :type include_sidechain: bool, optional
    :param peak_list: comma-separated two column file can be given as optional unassigned peak list,
        which can be overlaid on the spectrum as another data set, defaults to None
    :type peak_list: str, optional
    :param output_format: visualizations can be exported as interactive 'html' file
        or as static images in 'jpg','jpeg','png','pdf','webp','svg', defaults to 'html'
    :type output_format: str, optional
    :param output_file: file name to export visualization
    :type output_file: str, optional
    :param output_image_width: The width of the exported image in layout pixels, defaults to 800
    :type output_image_width: int, optional
    :param output_image_height: The height of the exported image in layout pixels, defaults to 600
    :type output_image_height: int, optional
    :param show_visualization: Visualization automatically opens in a web browser or as
        embedded visualization in Jupyter Notebook. This feature can be disabled
        by setting this flag as False, defaults to True
    :type show_visualization: bool, optional
    :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track)
        if draw_trace is True cs_track={ matching atoms:[cs_values]} else cs_track={}
    :rtype: tuple
    """

    if bmrb_ids is None and input_file_names is None and entry_objects is None:
        logging.error('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
        raise TypeError('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
    x1 = []
    y1 = []
    if peak_list is not None:
        with open(peak_list) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                x1.append(float(row[0]))
                y1.append(float(row[1]))
    peak_list_2d = create_n15hsqc_peaklist(bmrb_ids,
                                           input_file_names=input_file_names,
                                           entry_objects=entry_objects,
                                           auth_tag=auth_tag,
                                           draw_trace=draw_trace,
                                           include_sidechain=include_sidechain)
    x = peak_list_2d[0]
    y = peak_list_2d[1]
    data_set = peak_list_2d[2]
    info = peak_list_2d[3]
    res = peak_list_2d[4]
    cs_track = peak_list_2d[5]
    if (len(x) == 0 or len(y) == 0) and (len(x1) == 0 or len(y1) == 0):
        logging.error('Required chemical shifts not found')
        raise ValueError('Required chemical shifts not found')
    if legend is None:
        if len(x1) and len(y1):
            fig = px.scatter(x=x, y=y,
                             title='Simulated <sup>1</sup>H-<sup>15</sup>N HSQC peak positions',
                             hover_name=info,
                             color=data_set,
                             labels={"color": "Data set",
                                     # "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>15</sup>N (ppm)',
                                     }, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        else:
            fig = px.scatter(x=x, y=y,
                             title='Simulated <sup>1</sup>H-<sup>15</sup>N HSQC peak positions',
                             symbol=data_set,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>15</sup>N (ppm)',
                                     }, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        if len(x1) and len(y1):
            fig.update_layout(showlegend=True)
        else:
            fig.update_layout(showlegend=False)

        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")

    elif legend == 'residue':
        fig = px.scatter(x=x, y=y,
                         title='Simulated <sup>1</sup>H-<sup>15</sup>N HSQC peak positions',
                         hover_name=info,
                         color=res,
                         labels={"color": "Residue",
                                 # "symbol": "Data set",
                                 "x": '<sup>1</sup>H (ppm)',
                                 "y": '<sup>15</sup>N (ppm)',
                                 }, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, mode='lines',
                                opacity=0.7, hover_name=k)
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")

    elif legend == 'dataset':
        fig = px.scatter(x=x, y=y,
                         title='Simulated <sup>1</sup>H-<sup>15</sup>N HSQC peak positions',
                         hover_name=info,
                         color=data_set,
                         labels={"color": "Data set",
                                 # "symbol": "Data set",
                                 "x": '<sup>1</sup>H (ppm)',
                                 "y": '<sup>15</sup>N (ppm)',
                                 }, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0],
                                y=cs_track[k][1],
                                name=k,
                                opacity=0.7,
                                mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")
    else:
        raise ValueError('legend type not supported')

    if show_visualization:
        fig.show()
    if output_file is not None:
        if output_format == 'html':
            fig.write_html(output_file)
        elif output_format in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'json']:
            fig.write_image(file=output_file, format=output_format, width=output_image_width,
                            height=output_image_height)
        else:
            logging.error('Output file format not supported:{}'.format(output_format))
    return x, y, data_set, info, res, cs_track


def c13hsqc(bmrb_ids: Optional[Union[str, List[str], int, List[int]]] = None,
            input_file_names: Optional[Union[str, List[str]]] = None,
            entry_objects: Optional[Union[pynmrstar.Entry, List[pynmrstar.Entry]]] = None,
            auth_tag: Optional[bool] = False,
            legend: Optional[str] = None,
            draw_trace: Optional[bool] = False,
            peak_list: Optional[str] = None,
            output_format: Optional[str] = 'html',
            output_file: Optional[str] = None,
            output_image_width: Optional[int] = 800,
            output_image_height: Optional[int] = 600,
            show_visualization: Optional[bool] = True) -> tuple:
    """
    Plots \u00b9\u00b3C-HSQC spectrum  for a given  BMRB entry/NMR-STAR file/PyNMRSTAR entry object;
    This function can be used to compare different data sets as overlaid NMR Spectra. It overlays \u00b9\u00b3C-HSQC
    for a given list of BMRB entries/NMR-STAR files/PyNMRSTAR entry objects and draw lines connecting peaks
    from residues at the same sequence location in different data sets

    :param bmrb_ids: single BMRB entry ID or multiple BMRB entry IDs as list, defaults to None
    :type bmrb_ids: str/int/list, optional
    :param input_file_names: single input file name  or multiple input file names as list, defaults to None
    :type input_file_names: str, optional
    :param entry_objects: single PyNMRSTAR entry object  or multiple PyNMRSTAR entry objects as list, defaults to None
    :type entry_objects: PyNMRSTAR object/list, optional
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :param legend: legends are disabled by default. Residue types are indicated by  color and  data sets are
        indicated by symbol, displaying the combination of both will create a very long list of legend. Optionally
        either 'residue' or 'dataset' can be used to color code the scatter plot by residue type
        or data set and display the legend, defaults to None
    :type legend: str, optional
    :param draw_trace: draw line connecting peaks from residues at the same sequence location in different
        data sets, defaults to False
    :type draw_trace: bool optional
    :param peak_list: comma-separated two column file can be given as optional unassigned peak list,
        which can be overlaid on the spectrum as another data set, defaults to None
    :type peak_list: str, optional
    :param output_format: visualizations can be exported as interactive 'html' file
        or as static images in 'jpg','jpeg','png','pdf','webp','svg', defaults to 'html'
    :type output_format: str, optional
    :param output_file: file name to export visualization
    :type output_file: str, optional
    :param output_image_width: The width of the exported image in layout pixels, defaults to 800
    :type output_image_width: int, optional
    :param output_image_height: The height of the exported image in layout pixels, defaults to 600
    :type output_image_height: int, optional
    :param show_visualization: Visualization automatically opens in a web browser or as
        embedded visualization in Jupyter Notebook. This feature can be disabled
        by setting this flag as False, defaults to True
    :type show_visualization: bool, optional
    :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track)
        if draw_trace is True cs_track={ matching atoms:[cs_values]} else cs_track={}
    :rtype: tuple
    """
    if bmrb_ids is None and input_file_names is None and entry_objects is None:
        logging.error('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
        raise TypeError('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
    peak_list_2d = create_c13hsqc_peaklist(bmrb_ids,
                                           input_file_names=input_file_names,
                                           entry_objects=entry_objects,
                                           auth_tag=auth_tag,
                                           draw_trace=draw_trace)
    x1 = []
    y1 = []
    if peak_list is not None:
        with open(peak_list) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                x1.append(float(row[0]))
                y1.append(float(row[1]))

    x = peak_list_2d[0]
    y = peak_list_2d[1]
    data_set = peak_list_2d[2]
    info = peak_list_2d[3]
    res = peak_list_2d[4]
    cs_track = peak_list_2d[5]
    if len(x) == 0 or len(y) == 0:
        logging.error('Required chemical shifts not found')
        raise ValueError('Required chemical shifts not found')
    if legend is None:
        fig = px.scatter(x=x, y=y,
                         title='Simulated <sup>1</sup>H-<sup>13</sup>C HSQC peak positions',
                         symbol=data_set,
                         hover_name=info,
                         color=res,
                         labels={"color": "Residue",
                                 "symbol": "Data set",
                                 "x": '<sup>1</sup>H (ppm)',
                                 "y": '<sup>13</sup>C (ppm)'}, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1],
                                name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_layout(showlegend=False)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")

    elif legend == 'residue':
        fig = px.scatter(x=x, y=y,
                         title='Simulated <sup>1</sup>H-<sup>13</sup>C HSQC peak positions',
                         hover_name=info,
                         color=res,
                         labels={"color": "Residue",
                                 # "symbol": "Data set",
                                 "x": '<sup>1</sup>H (ppm)',
                                 "y": '<sup>13</sup>C (ppm)'}, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, mode='lines', opacity=0.7)
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")

    elif legend == 'dataset':
        fig = px.scatter(x=x, y=y,
                         title='Simulated <sup>1</sup>H-<sup>13</sup>C HSQC peak positions',
                         hover_name=info,
                         color=data_set,
                         labels={"color": "Data set",
                                 # "symbol": "Data set",
                                 "x": '<sup>1</sup>H (ppm)',
                                 "y": '<sup>13</sup>C (ppm)'}, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")
    else:
        raise ValueError('legend type not supported')

    if show_visualization:
        fig.show()
    if output_file is not None:
        if output_format == 'html':
            fig.write_html(output_file)
        elif output_format in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'json']:
            fig.write_image(file=output_file, format=output_format, width=output_image_width,
                            height=output_image_height)
        else:
            logging.error('Output file format not supported:{}'.format(output_format))
    return x, y, data_set, info, res, cs_track


def tocsy(bmrb_ids: Optional[Union[str, List[str], int, List[int]]] = None,
          input_file_names: Optional[Union[str, List[str]]] = None,
          entry_objects: Optional[Union[pynmrstar.Entry, List[pynmrstar.Entry]]] = None,
          auth_tag: Optional[bool] = False,
          legend: Optional[str] = None,
          draw_trace: Optional[bool] = False,
          peak_list: Optional[str] = None,
          output_format: Optional[str] = 'html',
          output_file: Optional[str] = None,
          output_image_width: Optional[int] = 800,
          output_image_height: Optional[int] = 600,
          show_visualization: Optional[bool] = True) -> tuple:
    """
    Plots \u00b9H-\u00b9H-TOCSY spectrum for a given  BMRB entry/NMR-STAR file/PyNMRSTAR entry object;
    This function can be used to compare different data sets as overlaid NMR Spectra. It overlays \u00b9H-\u00b9H-TOCSY
    for a given list of BMRB entries/NMR-STAR files/PyNMRSTAR entry objects and draw lines connecting peaks
    from residues at the same sequence location in different data sets

    :param bmrb_ids: single BMRB entry ID or multiple BMRB entry IDs as list, defaults to None
    :type bmrb_ids: str/int/list, optional
    :param input_file_names: single input file name  or multiple input file names as list, defaults to None
    :type input_file_names: str, optional
    :param entry_objects: single PyNMRSTAR entry object  or multiple PyNMRSTAR entry objects as list, defaults to None
    :type entry_objects: PyNMRSTAR object/list, optional
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :param legend: legends are disabled by default. Residue types are indicated by  color and  data sets are
        indicated by symbol, displaying the combination of both will create a very long list of legend. Optionally
        either 'residue' or 'dataset' can be used to color code the scatter plot by residue type
        or data set and display the legend, defaults to None
    :type legend: str, optional
    :param draw_trace: draw line connecting peaks from residues at the same sequence location in different
        data sets, defaults to False
    :type draw_trace: bool optional
    :param peak_list: comma-separated two column file can be given as optional unassigned peak list,
        which can be overlaid on the spectrum as another data set, defaults to None
    :type peak_list: str, optional
    :param output_format: visualizations can be exported as interactive 'html' file
        or as static images in 'jpg','jpeg','png','pdf','webp','svg', defaults to 'html'
    :type output_format: str, optional
    :param output_file: file name to export visualization
    :type output_file: str, optional
    :param output_image_width: The width of the exported image in layout pixels, defaults to 800
    :type output_image_width: int, optional
    :param output_image_height: The height of the exported image in layout pixels, defaults to 600
    :type output_image_height: int, optional
    :param show_visualization: Visualization automatically opens in a web browser or as
        embedded visualization in Jupyter Notebook. This feature can be disabled
        by setting this flag as False, defaults to True
    :type show_visualization: bool, optional
    :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track)
        if draw_trace is True cs_track={ matching atoms:[cs_values]} else cs_track={}
    :rtype: tuple
    """
    if bmrb_ids is None and input_file_names is None and entry_objects is None:
        logging.error('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
        raise TypeError('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')

    peak_list_2d = create_tocsy_peaklist(bmrb_ids,
                                         input_file_names=input_file_names,
                                         entry_objects=entry_objects,
                                         auth_tag=auth_tag,
                                         draw_trace=draw_trace)
    x1 = []
    y1 = []
    if peak_list is not None:
        with open(peak_list) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                x1.append(float(row[0]))
                y1.append(float(row[1]))
    x = peak_list_2d[0]
    y = peak_list_2d[1]
    data_set = peak_list_2d[2]
    info = peak_list_2d[3]
    res = peak_list_2d[4]
    cs_track = peak_list_2d[5]
    if len(x) == 0 or len(y) == 0:
        logging.error('Required chemical shifts not found')
        raise ValueError('Required chemical shifts not found')
    if legend is None:
        fig = px.scatter(x=x, y=y,
                         title='Simulated <sup>1</sup>H-<sup>1</sup>H TOCSY peak positions',
                         symbol=data_set,
                         hover_name=info,
                         color=res,
                         labels={"color": "Residue",
                                 "symbol": "Data set",
                                 "x": '<sup>1</sup>H (ppm)',
                                 "y": '<sup>1</sup>H (ppm)'}, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_layout(showlegend=False)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")

    elif legend == 'residue':
        fig = px.scatter(x=x, y=y,
                         title='Simulated <sup>1</sup>H-<sup>1</sup>H TOCSY peak positions',
                         hover_name=info,
                         color=res,
                         labels={"color": "Residue",
                                 # "symbol": "Data set",
                                 "x": '<sup>1</sup>H (ppm)',
                                 "y": '<sup>1</sup>H (ppm)'}, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, mode='lines',
                                opacity=0.7)
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")

    elif legend == 'dataset':
        fig = px.scatter(x=x, y=y,
                         title='Simulated <sup>1</sup>H-<sup>1</sup>H TOCSY peak positions',
                         hover_name=info,
                         color=data_set,
                         labels={"color": "Data set",
                                 # "symbol": "Data set",
                                 "x": '<sup>1</sup>H (ppm)',
                                 "y": '<sup>1</sup>H (ppm)'}, opacity=0.7).update(layout=dict(title=dict(x=0.5)))
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")
    else:
        raise ValueError('legend type not supported')
    if show_visualization:
        fig.show()
    if output_file is not None:
        if output_format == 'html':
            fig.write_html(output_file)
        elif output_format in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'json']:
            fig.write_image(file=output_file, format=output_format, width=output_image_width,
                            height=output_image_height)
        else:
            logging.error('Output file format not supported:{}'.format(output_format))
    return x, y, data_set, info, res, cs_track


def generic_2d(atom_x: str,
               atom_y: str,
               bmrb_ids: Optional[Union[str, List[str], int, List[int]]] = None,
               input_file_names: Optional[Union[str, List[str]]] = None,
               entry_objects: Optional[Union[pynmrstar.Entry, List[pynmrstar.Entry]]] = None,
               auth_tag: Optional[bool] = False,
               legend: Optional[str] = None,
               draw_trace: Optional[bool] = False,
               peak_list: Optional[str] = None,
               include_preceding: Optional[bool] = False,
               include_next: Optional[bool] = False,
               full_walk: Optional[bool] = False,
               seq_walk: Optional[bool] = False,
               output_format: Optional[str] = 'html',
               output_file: Optional[str] = None,
               output_image_width: Optional[int] = 800,
               output_image_height: Optional[int] = 600,
               show_visualization: Optional[bool] = True) -> tuple:
    """
    Plots generic 2D spectrum for a given  BMRB entry/NMR-STAR file/PyNMRSTAR entry object;
    This function can be used to compare different data sets as overlaid NMR Spectra. It overlays the 2D spectra
    for a given list of BMRB entries/NMR-STAR files/PyNMRSTAR entry objects and draw lines connecting peaks
    from residues at the same sequence location in different data sets

    :param atom_x: atom name for X coordinate in IUPAC format
    :type atom_x: str
    :param atom_y: atom name for Y coordinate in IUPAC format
    :type atom_y: str
    :param bmrb_ids: single BMRB entry ID or multiple BMRB entry IDs as list, defaults to None
    :type bmrb_ids: str/int/list, optional
    :param input_file_names: single input file name  or multiple input file names as list, defaults to None
    :type input_file_names: str, optional
    :param entry_objects: single PyNMRSTAR entry object  or multiple PyNMRSTAR entry objects as list, defaults to None
    :type entry_objects: PyNMRSTAR object/list, optional
    :param auth_tag: Use sequence numbers from _Atom_chem_shift.Auth_seq_ID instead of _Atom_chem_shift.Comp_index_ID
        in the NMR-STAR file/BMRB entry, defaults to False
    :type auth_tag: bool, optional
    :param legend: legends are disabled by default. Residue types are indicated by  color and  data sets are
        indicated by symbol, displaying the combination of both will create a very long list of legend. Optionally
        either 'residue' or 'dataset' can be used to color code the scatter plot by residue type
        or data set and display the legend, defaults to None
    :type legend: str, optional
    :param draw_trace: draw line connecting peaks from residues at the same sequence location in different
        data sets, defaults to False
    :type draw_trace: bool optional
    :param include_preceding: include peaks from i-1 residue on the Y axis, defaults to False
    :type include_preceding: bool, optional
    :param include_next: include peaks from i+i residue on the Y axis, defaults to False
    :type include_next: bool, optional
    :param seq_walk: draw line connecting i->i-1/i+1 to next i->i-1/i+1 for only
        continuous sequence segments, defaults to False
    :type seq_walk: bool, optional
    :param full_walk: draw line connecting i->i-1/i+1 to next i->i-1/i+1 for the
        full sequence ignoring any missing residues, defaults to False
    :type full_walk: bool, optional
    :param peak_list: comma-separated two column file can be given as optional unassigned peak list,
        which can be overlaid on the spectrum as another data set, defaults to None
    :type peak_list: str, optional
    :param output_format: visualizations can be exported as interactive 'html' file
        or as static images in 'jpg','jpeg','png','pdf','webp','svg', defaults to 'html'
    :type output_format: str, optional
    :param output_file: file name to export visualization
    :type output_file: str, optional
    :param output_image_width: The width of the exported image in layout pixels, defaults to 800
    :type output_image_width: int, optional
    :param output_image_height: The height of the exported image in layout pixels, defaults to 600
    :type output_image_height: int, optional
    :param show_visualization: Visualization automatically opens in a web browser or as
        embedded visualization in Jupyter Notebook. This feature can be disabled
        by setting this flag as False, defaults to True
    :type show_visualization: bool, optional
    :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track)
        if draw_trace is True cs_track={ matching atoms:[cs_values]} else cs_track={}
    :rtype: tuple
    """
    if bmrb_ids is None and input_file_names is None and entry_objects is None:
        logging.error('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
        raise TypeError('At least one of  three parameters must be present; bmrb_ids, input_file_names, entry_objects')
    peak_list_2d = create_2d_peaklist(atom_x=atom_x,
                                      atom_y=atom_y,
                                      bmrb_ids=bmrb_ids,
                                      input_file_names=input_file_names,
                                      entry_objects=entry_objects,
                                      auth_tag=auth_tag,
                                      draw_trace=draw_trace,
                                      include_preceding=include_preceding,
                                      include_next=include_next,
                                      legend=legend)
    x1 = []
    y1 = []
    if peak_list is not None:
        with open(peak_list) as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                x1.append(float(row[0]))
                y1.append(float(row[1]))
    x = peak_list_2d[0]
    y = peak_list_2d[1]
    data_set = peak_list_2d[2]
    info = peak_list_2d[3]
    res = peak_list_2d[4]
    cs_track = peak_list_2d[5]
    psn = peak_list_2d[6]
    seq_trace = peak_list_2d[7]

    def sequence_walk(seq_t, filt=True):
        seq_traces = {}
        j = 0
        s = -999
        for i in range(len(seq_t[0])):
            if seq_t[2][i] == s or seq_t[2][i] == s + 1:
                s = seq_t[2][i]
                seq_traces['walk_{}'.format(j)][0].append(seq_t[0][i])
                seq_traces['walk_{}'.format(j)][1].append(seq_t[1][i])
            else:
                s = seq_t[2][i]
                j += 1
                seq_traces['walk_{}'.format(j)] = [[seq_t[0][i]], [seq_t[1][i]]]
        sq_w = {}
        if filt:
            for k1 in seq_traces.keys():
                if len(seq_traces[k1][0]) > 1:
                    sq_w[k1] = seq_traces[k1]
        else:
            sq_w = seq_traces
        return sq_w

    sq_walk = sequence_walk(seq_trace)
    if len(x) == 0 or len(y) == 0:
        logging.error('Required chemical shifts not found')
        raise ValueError('Required chemical shifts not found')
    if legend is None:
        fig = px.scatter(x=x, y=y,
                         title='Simulated {}-{} COSY peak positions'.format(atom_x, atom_y),
                         symbol=data_set,
                         hover_name=info,
                         size=psn,
                         color=res,
                         labels={"color": "Residue",
                                 "symbol": "Data set",
                                 "x": '{} (ppm)'.format(atom_x),
                                 "y": '{} (ppm)'.format(atom_y)}, opacity=0.7)
        fig.update(layout=dict(title=dict(x=0.5)))
        if full_walk:
            fig.add_scatter(x=seq_trace[0], y=seq_trace[1], name='Full walk', opacity=0.7, mode='lines')
        if seq_walk:
            for k in sq_walk.keys():
                fig.add_scatter(x=sq_walk[k][0], y=sq_walk[k][1], name=k, opacity=0.7, mode='lines')
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_layout(showlegend=False)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")

    elif legend == 'residue':
        fig = px.scatter(x=x, y=y,
                         title='Simulated {}-{} COSY peak positions'.format(atom_x, atom_y),
                         hover_name=info,
                         color=res,
                         size=psn,
                         labels={"color": "Residue",
                                 # "symbol": "Data set",
                                 "x": '{} (ppm)'.format(atom_x),
                                 "y": '{} (ppm)'.format(atom_y)}, opacity=0.7).update(
            layout=dict(title=dict(x=0.5)))
        if full_walk:
            fig.add_scatter(x=seq_trace[0], y=seq_trace[1], name='Full walk', opacity=0.7, mode='lines')
        if seq_walk:
            for k in sq_walk.keys():
                fig.add_scatter(x=sq_walk[k][0], y=sq_walk[k][1], name=k, opacity=0.7, mode='lines')
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")

    elif legend == 'dataset':
        fig = px.scatter(x=x, y=y,
                         title='Simulated {}-{} COSY peak positions'.format(atom_x, atom_y),
                         hover_name=info,
                         color=data_set,
                         size=psn,
                         labels={"color": "Data set",
                                 # "symbol": "Data set",
                                 "x": '{} (ppm)'.format(atom_x),
                                 "y": '{} (ppm)'.format(atom_y)}, opacity=0.7).update(
            layout=dict(title=dict(x=0.5)))
        if full_walk:
            fig.add_scatter(x=seq_trace[0], y=seq_trace[1], name='Full walk', opacity=0.7, mode='lines')
        if seq_walk:
            for k in sq_walk.keys():
                fig.add_scatter(x=sq_walk[k][0], y=sq_walk[k][1], name=k, opacity=0.7, mode='lines')
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")
    elif legend == 'psn':
        fig = px.scatter(x=x, y=y,
                         title='Simulated {}-{} COSY peak positions'.format(atom_x, atom_y),
                         hover_name=info,
                         color=psn,
                         labels={"color": "Data set",
                                 # "symbol": "Data set",
                                 "x": '{} (ppm)'.format(atom_x),
                                 "y": '{} (ppm)'.format(atom_y)}, opacity=0.7).update(
            layout=dict(title=dict(x=0.5)))
        if full_walk:
            fig.add_scatter(x=seq_trace[0], y=seq_trace[1], name='Full walk', opacity=0.7, mode='lines')
        if seq_walk:
            for k in sq_walk.keys():
                fig.add_scatter(x=sq_walk[k][0], y=sq_walk[k][1], name=k, opacity=0.7, mode='lines')
        if draw_trace:
            for k in cs_track.keys():
                fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7, mode='lines')
        if peak_list is not None:
            fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")
    else:
        raise ValueError('legend type not supported')
    if show_visualization:
        fig.show()
    if output_file is not None:
        if output_format == 'html':
            fig.write_html(output_file)
        elif output_format in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'json']:
            fig.write_image(file=output_file, format=output_format, width=output_image_width,
                            height=output_image_height)
        else:
            logging.error('Output file format not supported:{}'.format(output_format))
    return x, y, data_set, info, res, cs_track


def export_peak_list(peak_list: tuple,
                     output_file_name: Optional[str] = None,
                     output_format: str = 'csv',
                     include_side_chain: bool = True) -> dict:
    """
    Exports the peak list return from other functions in this module in csv or sparky format

    :param peak_list: Output tuple from any peak list/spectrum simulation function from the module Spectra
    :type peak_list: tuple
    :param output_file_name: output file name
    :type output_file_name: str, optional
    :param output_format: output format 'csv' or 'sparky', defaults to 'csv'
    :type output_format: str, optional
    :param include_side_chain: whether or not include side chain resonances in the output, defaults to True
    :type include_side_chain: bool, optional
    :return: data dictionary {'column header1':[values],'column header1':[values]..}
    :rtype: dict
    """
    back_bone = ['H', 'N', 'C', 'CA']
    if output_format == 'csv':
        csv_dict = {'sequence': [],
                    'chem_comp_ID': [],
                    'X_shift': [],
                    'Y_shift': [],
                    'X_atom_name': [],
                    'Y_atom_name': []}
        for i in range(len(peak_list[0])):
            atom_x = peak_list[3][i].split("-")[5]
            atom_y = peak_list[3][i].split("-")[6]
            if not include_side_chain:
                if atom_x in back_bone and atom_y in back_bone:
                    csv_dict['sequence'].append(peak_list[3][i].split("-")[3])
                    csv_dict['chem_comp_ID'].append(peak_list[3][i].split("-")[4])
                    csv_dict['X_shift'].append(peak_list[0][i])
                    csv_dict['Y_shift'].append(peak_list[1][i])
                    csv_dict['X_atom_name'].append(peak_list[3][i].split("-")[5])
                    csv_dict['Y_atom_name'].append(peak_list[3][i].split("-")[6])

            else:
                csv_dict['sequence'].append(peak_list[3][i].split("-")[3])
                csv_dict['chem_comp_ID'].append(peak_list[3][i].split("-")[4])
                csv_dict['X_shift'].append(peak_list[0][i])
                csv_dict['Y_shift'].append(peak_list[1][i])
                csv_dict['X_atom_name'].append(peak_list[3][i].split("-")[5])
                csv_dict['Y_atom_name'].append(peak_list[3][i].split("-")[6])
        if output_file_name is not None:
            fo = open(output_file_name, 'w')
            fo.write('sequence,chem_comp_ID,X_shift,Y_shift,X_atom_name,Y_atom_name\n')
            for i in range(len(csv_dict['sequence'])):
                fo.write('{},{},{},{},{},{}\n'.format(csv_dict['sequence'][i],
                                                      csv_dict['chem_comp_ID'][i],
                                                      round(float(csv_dict['X_shift'][i]), 3),
                                                      round(float(csv_dict['Y_shift'][i]), 3),
                                                      csv_dict['X_atom_name'][i],
                                                      csv_dict['Y_atom_name'][i]))
            fo.close()
    elif output_format == 'sparky':
        csv_dict = {'Assignment': [],
                    'w1': [],
                    'w2': []}
        for i in range(len(peak_list[0])):
            atom_x = peak_list[3][i].split("-")[5]
            atom_y = peak_list[3][i].split("-")[6]
            res = peak_list[3][i].split("-")[4]
            try:
                assignment = '{}{}{}-{}'.format(ChemicalShiftStatistics.one_letter_code[res],
                                                peak_list[3][i].split("-")[3],
                                                atom_x, atom_y)
            except KeyError:
                assignment = '{}{}{}-{}'.format('X',
                                                peak_list[3][i].split("-")[3],
                                                atom_x, atom_y)
            if not include_side_chain:
                if atom_x in back_bone and atom_y in back_bone:
                    csv_dict['Assignment'].append(assignment)
                    csv_dict['w1'].append(peak_list[1][i])
                    csv_dict['w2'].append(peak_list[0][i])
            else:
                csv_dict['Assignment'].append(assignment)
                csv_dict['w1'].append(peak_list[1][i])
                csv_dict['w2'].append(peak_list[0][i])
        if output_file_name is not None:
            fo = open(output_file_name, 'w')
            fo.write('Assignment  \t{:>6}\t\t{:>6}\n\n'.format('w1', 'w2'))
            for i in range(len(csv_dict['Assignment'])):
                fo.write('{:<10}\t\t{:>6}\t\t{:>6}\n'.format(csv_dict['Assignment'][i],
                                                             round(float(csv_dict['w1'][i]), 3),
                                                             round(float(csv_dict['w2'][i]), 3)))
            fo.close()
    else:
        logging.error('Output format not supported')
        raise ValueError('Output format not supported')
    return csv_dict

#
# if __name__ == "__main__":
#
#     # The following script is to generate example figures for readthedocs
#     n15hsqc(bmrb_ids=[17074,17076,17077],legend='dataset',output_format='jpg',
#             output_file='../docs/_images/sample_n15hsqc.jpg',show_visualization=False,
#             draw_trace=True)
#     n15hsqc(bmrb_ids=[17074, 17076, 17077], legend='dataset', output_format='html',
#             output_file='../docs/_static/sample_n15hsqc.html', show_visualization=False,
#             draw_trace=True)
#     c13hsqc(bmrb_ids=15060,legend='residue',output_file='../docs/_images/sample_c13hsqc.jpg',
#             output_format='jpg',show_visualization=False)
#     c13hsqc(bmrb_ids=15060, legend='residue', output_file='../docs/_static/sample_c13hsqc.html',
#             output_format='html', show_visualization=False)
#     n15hsqc(bmrb_ids=[17076, 17077], input_file_names='tests/test_data/MyData.str',
#                                 legend='dataset',output_format='jpg',
#             output_file='../docs/_images/quick_start_n15hsqc_compare.jpg',
#             show_visualization=False)
#     n15hsqc(bmrb_ids=[17076, 17077], input_file_names='tests/test_data/MyData.str',
#             legend='dataset', output_format='html',
#             output_file='../docs/_static/quick_start_n15hsqc_compare.html',
#             show_visualization=False)
#     n15hsqc(bmrb_ids=[17076, 17077], input_file_names='tests/test_data/MyData.str',
#             legend='dataset', output_format='html',draw_trace=True,
#             output_file='../docs/_static/quick_start_n15hsqc_compare2.html',
#             show_visualization=False)
#     n15hsqc(bmrb_ids=[17076, 17077], peak_list='tests/test_data/test_peak_list.csv',
#                                 legend='dataset', output_format='jpg',
#             output_file='../docs/_images/quick_start_n15_peaklist.jpg',
#             show_visualization=False)
#     n15hsqc(bmrb_ids=15060,
#             legend='residue',
#             output_format='jpg',
#             output_file='../docs/_images/example1.jpg',
#             show_visualization=False)
#     n15hsqc(bmrb_ids=15060,
#             legend='residue',
#             output_format='html',
#             output_file='../docs/_static/example1.html',
#             show_visualization=False)
#     n15hsqc(bmrb_ids=[17076, 17077],
#             input_file_names='tests/test_data/MyData.str',
#             legend='dataset',
#             output_format='jpg',
#             output_file='../docs/_images/example2.jpg',
#             show_visualization=False
#             )
#     n15hsqc(bmrb_ids=[17076, 17077],
#             input_file_names='tests/test_data/MyData.str',
#             legend='dataset',
#             output_format='html',
#             output_file='../docs/_static/example2.html',
#             show_visualization=False
#             )
#     n15hsqc(bmrb_ids=[17076, 17077],
#             input_file_names='tests/test_data/MyData.str',
#             peak_list='tests/test_data/test_peak_list.csv',
#             legend='dataset',
#             output_format='jpg',
#             output_file='../docs/_images/example3.jpg',
#             show_visualization=False
#             )
#     n15hsqc(bmrb_ids=[17076, 17077],
#             input_file_names='tests/test_data/MyData.str',
#             peak_list='tests/test_data/test_peak_list.csv',
#             legend='dataset',
#             output_format='html',
#             output_file='../docs/_static/example3.html',
#             show_visualization=False
#             )
#     n15hsqc(bmrb_ids=[17076, 17077],
#             input_file_names='tests/test_data/MyData.str',
#             legend='dataset',
#             draw_trace=True,
#             output_format='jpg',
#             output_file='../docs/_images/example4.jpg',
#             show_visualization=False
#             )
#     n15hsqc(bmrb_ids=[17076, 17077],
#             input_file_names='tests/test_data/MyData.str',
#             legend='dataset',
#             draw_trace=True,
#             output_format='html',
#             output_file='../docs/_static/example4.html',
#             show_visualization=False
#             )
#     c13hsqc(bmrb_ids=15060,
#             legend='residue',
#             output_format='jpg',
#             output_file='../docs/_images/example5.jpg',
#             show_visualization=False)
#     c13hsqc(bmrb_ids=15060,
#             legend='residue',
#             output_format='html',
#             output_file='../docs/_static/example5.html',
#             show_visualization=False)
#     c13hsqc(bmrb_ids=[17074, 17076, 17077],
#             legend='dataset',
#             output_format='jpg',
#             output_file='../docs/_images/example6.jpg',
#             show_visualization=False
#             )
#     c13hsqc(bmrb_ids=[17074, 17076, 17077],
#             legend='dataset',
#             output_format='html',
#             output_file='../docs/_static/example6.html',
#             show_visualization=False
#             )
#     c13hsqc(bmrb_ids=[17074, 17076, 17077],
#             legend='dataset',
#             draw_trace=True,
#             output_format='jpg',
#             output_file='../docs/_images/example7.jpg',
#             show_visualization=False
#             )
#     c13hsqc(bmrb_ids=[17074, 17076, 17077],
#             legend='dataset',
#             draw_trace=True,
#             output_format='html',
#             output_file='../docs/_static/example7.html',
#             show_visualization=False
#             )
#     tocsy(bmrb_ids=15060,
#           legend='residue',
#           output_format='jpg',
#           output_file='../docs/_images/example8.jpg',
#           show_visualization=False)
#     tocsy(bmrb_ids=15060,
#           legend='residue',
#           output_format='html',
#           output_file='../docs/_static/example8.html',
#           show_visualization=False)
#     tocsy(bmrb_ids=[17074, 17076, 17077],
#           legend='dataset',
#           output_format='jpg',
#           output_file='../docs/_images/example9.jpg',
#           show_visualization=False)
#     tocsy(bmrb_ids=[17074, 17076, 17077],
#           legend='dataset',
#           output_format='html',
#           output_file='../docs/_static/example9.html',
#           show_visualization=False)
#     tocsy(bmrb_ids=[17074, 17076, 17077],
#           legend='residue',
#           output_format='jpg',
#           output_file='../docs/_images/example10.jpg',
#           show_visualization=False)
#     tocsy(bmrb_ids=[17074, 17076, 17077],
#           legend='residue',
#           output_format='html',
#           output_file='../docs/_static/example10.html',
#           show_visualization=False)
#     tocsy(bmrb_ids=[17074, 17076, 17077],
#           legend='dataset',
#           draw_trace=True,
#           output_format='jpg',
#           output_file='../docs/_images/example11.jpg',
#           show_visualization=False)
#     tocsy(bmrb_ids=[17074, 17076, 17077],
#           legend='dataset',
#           draw_trace=True,
#           output_format='html',
#           output_file='../docs/_static/example11.html',
#           show_visualization=False)
#     generic_2d(bmrb_ids=15060,
#                atom_x='N',
#                atom_y='CB',
#                legend='residue',
#                output_format='jpg',
#                output_file='../docs/_images/example12.jpg',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15060,
#                atom_x='N',
#                atom_y='CB',
#                legend='residue',
#                output_format='html',
#                output_file='../docs/_static/example12.html',
#                show_visualization=False)
#     generic_2d(bmrb_ids=[17074,17076,17077],
#                atom_x='N',
#                atom_y='CB',
#                legend='dataset',
#                output_format='jpg',
#                output_file='../docs/_images/example13.jpg',
#                show_visualization=False)
#     generic_2d(bmrb_ids=[17074,17076,17077],
#                atom_x='N',
#                atom_y='CB',
#                legend='dataset',
#                output_format='html',
#                output_file='../docs/_static/example13.html',
#                show_visualization=False)
#     generic_2d(bmrb_ids=[17074, 17076, 17077],
#                atom_x='N',
#                atom_y='CB',
#                legend='dataset',
#                draw_trace=True,
#                output_format='jpg',
#                output_file='../docs/_images/example14.jpg',
#                show_visualization=False)
#     generic_2d(bmrb_ids=[17074, 17076, 17077],
#                atom_x='N',
#                atom_y='CB',
#                legend='dataset',
#                draw_trace=True,
#                output_format='html',
#                output_file='../docs/_static/example14.html',
#                show_visualization=False)
#
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_preceding=True,
#                output_format='jpg',
#                output_file='../docs/_images/example15.jpg',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_preceding=True,
#                output_format='html',
#                output_file='../docs/_static/example15.html',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_preceding=True,
#                seq_walk=True,
#                output_format='jpg',
#                output_file='../docs/_images/example16.jpg',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_preceding=True,
#                seq_walk=True,
#                output_format='html',
#                output_file='../docs/_static/example16.html',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_next=True,
#                output_format='jpg',
#                output_file='../docs/_images/example17.jpg',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_next=True,
#                output_format='html',
#                output_file='../docs/_static/example17.html',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_next=True,
#                seq_walk=True,
#                output_format='jpg',
#                output_file='../docs/_images/example18.jpg',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_next=True,
#                seq_walk=True,
#                output_format='html',
#                output_file='../docs/_static/example18.html',
#                show_visualization=False)
#
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_next=True,
#                full_walk=True,
#                output_format='jpg',
#                output_file='../docs/_images/example19.jpg',
#                show_visualization=False)
#     generic_2d(bmrb_ids=15000,
#                atom_x='N',
#                atom_y='CA',
#                legend='residue',
#                include_next=True,
#                full_walk=True,
#                output_format='html',
#                output_file='../docs/_static/example19.html',
#                show_visualization=False)
#

from __future__ import print_function

import json
import logging
import ntpath
import optparse
import os
import sys
import re
import numpy as np
import plotly
import pynmrstar

# Set the log level to INFO
logging.getLogger().setLevel(logging.INFO)

# Determine if we are running in python3
PY3 = (sys.version_info[0] == 3)
(scriptPath, scriptName) = ntpath.split(os.path.realpath(__file__))
# pylint: disable=wrong-import-position,no-name-in-module
# pylint: disable=import-error,wrong-import-order
# Python version dependent loads
if PY3:
    from urllib.request import urlopen, Request
else:
    from urllib2 import urlopen, Request

# _API_URL = "http://api.bmrb.io/v2"
_API_URL = "http://dev-api.bmrb.io/v2"
NOTEBOOK = False
_OPACITY = 0.5
_AUTOOPEN = True
__version__ = "test"

three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                     'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                     'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                     'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])


class ChemicalShift(object):

    def __init__(self):
        return

    @staticmethod
    def from_entry(entry_data, data_set_id, auth_tag=False, ):
        '''
        Extracts chemical shift data as dictionary from pynmrstar entry object
        :param entry_data: PyNMRSTAR entry object
        :param data_set_id: Data set identifier (bmrb id or filename or user defined id)
        :param auth_tag: Use author sequence numbering True/False default: False
        :return: Chemical shift dictionary
        '''
        cs_data = {}
        # entry_data = pynmrstar.Entry.from_database(bmrb_id)
        cs_loops = entry_data.get_loops_by_category('Atom_chem_shift')
        cs_list_id = 1
        for cs_loop in cs_loops:
            cs_id = '{}-{}'.format(data_set_id, cs_list_id)
            cs_list_id += 1
            cs_data[cs_id] = {}
            column_name = cs_loop.get_tag_names()
            data = cs_loop.data
            if auth_tag:
                # chain = column_name.index('_Atom_chem_shift.Auth_asym_ID')
                seq_no = column_name.index('_Atom_chem_shift.Auth_seq_ID')
                # residue = column_name.index('_Atom_chem_shift.Auth_comp_ID')
                # atom = column_name.index('_Atom_chem_shift.Auth_atom_ID')
            else:
                # chain = column_name.index('_Atom_chem_shift.Entity_assembly_ID')
                seq_no = column_name.index('_Atom_chem_shift.Comp_index_ID')
                # residue = column_name.index('_Atom_chem_shift.Comp_ID')
                # atom = column_name.index('_Atom_chem_shift.Atom_ID')
            residue = column_name.index('_Atom_chem_shift.Comp_ID')
            atom = column_name.index('_Atom_chem_shift.Atom_ID')
            chain = column_name.index('_Atom_chem_shift.Entity_assembly_ID')
            cs_value = column_name.index('_Atom_chem_shift.Val')
            for cs_row in data:
                c = cs_row[chain]
                s = int(cs_row[seq_no])
                r = cs_row[residue]
                a = cs_row[atom]
                v = float(cs_row[cs_value])
                if c not in cs_data[cs_id].keys():
                    cs_data[cs_id][c] = {}
                if s not in cs_data[cs_id][c].keys():
                    cs_data[cs_id][c][s] = {}
                if a not in cs_data[cs_id][c][s].keys():
                    cs_data[cs_id][c][s][a] = []
                cs_data[cs_id][c][s][a].append(r)
                try:
                    cs_data[cs_id][c][s][a].append(one_letter_code[r])
                except KeyError:
                    cs_data[cs_id][c][s][a].append(r)
                # cs_data[cs_id][c][s][a].append(a)
                cs_data[cs_id][c][s][a].append(v)
        for cs_id in cs_data.keys():
            for chain in cs_data[cs_id].keys():
                seq_ids = sorted(set([i for i in cs_data[cs_id][chain].keys()]))
                cs_data[cs_id][chain]['seq_ids'] = seq_ids
        return cs_data

    @classmethod
    def from_file(self, file_names, auth_tag=False, data_set_id=None):
        '''
        Extracts chemical shift information from NMR-STAR file
        :param file_name: list of NMR-STAR file names with full path or single file name
        :param auth_tag: Use author sequence numbering True/False default: False
        :param data_set_id: User defined data set id default: filename
        :return: Chemical shift dictionary
        '''
        # handle the data_set_id list for more than one files : todo
        if type(file_names) is list:
            all_cs_data = {}
            for file_name in file_names:
                if os.path.exists(file_name):
                    entry_data = pynmrstar.Entry.from_file(file_name)
                    data_set_id = os.path.splitext(os.path.basename(file_name))[0]
                    cs_data = self.from_entry(entry_data, data_set_id, auth_tag)
                else:
                    raise IOError('File not found : {}'.format(file_name))
                all_cs_data.update(cs_data)
        else:
            if os.path.exists(file_names):
                entry_data = pynmrstar.Entry.from_file(file_names)
                if data_set_id is None:
                    data_set_id = os.path.splitext(os.path.basename(file_names))[0]
                all_cs_data = self.from_entry(entry_data, data_set_id, auth_tag)
            else:
                raise IOError('File not found : {}'.format(file_names))
        return all_cs_data

    @classmethod
    def from_bmrb(self, bmrb_ids, auth_tag=False):
        '''
        Extracts chemical shift information directly from BMRB database for a given BMRB entry id
        :param bmrb_ids: list of BMRB entry ids or single BMRB ID
        :param auth_tag: Use author sequence numbering True/False default: False
        :return: Chemical shift dictionary
        '''
        if type(bmrb_ids) is list:
            all_cs_data={}
            for bmrb_id in bmrb_ids:
                try:
                    entry_data = pynmrstar.Entry.from_database(bmrb_id)
                except OSError:
                    entry_data = None
                except KeyError:
                    entry_data = None
                except IOError:
                    entry_data = None
                if entry_data is not None:
                    cs_data = self.from_entry(entry_data, bmrb_id, auth_tag)
                else:
                    raise IOError('Entry not found in public database: {}'.format(bmrb_id))
                all_cs_data.update(cs_data)
        else:
            try:
                entry_data = pynmrstar.Entry.from_database(bmrb_ids)
            except OSError:
                entry_data = None
            except KeyError:
                entry_data = None
            except IOError:
                entry_data = None
            if entry_data is not None:
                all_cs_data = self.from_entry(entry_data, bmrb_ids, auth_tag)
            else:
                raise IOError('Entry not found in public database: {}'.format(bmrb_id))

        return all_cs_data


    def from_bmrb_database(self):
        pass


if __name__ == "__main__":
    p = ChemicalShift()
    p.from_file('/Users/kumaran/MyData3.str',data_set_id='test')

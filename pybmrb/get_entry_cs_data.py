#!/usr/bin/env python3

from __future__ import print_function
import logging
import os
import pynmrstar
from typing import TextIO, BinaryIO, Union, List, Optional, Dict, Any, Tuple

# Set the log level to INFO
logging.getLogger().setLevel(logging.INFO)

three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                     'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                     'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                     'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])


class ChemicalShift(object):
    """
    Generates chemical shift dictionary from BMRB entry or NMR-STAR file
    """

    @staticmethod
    def _from_pynmrstar_entry_object(entry_data: pynmrstar.Entry,
                                     data_set_id: str,
                                     auth_tag: bool = False, ):
        """
        Extracts chemical shift data as dictionary from PyNMRSTAR entry object

        :param entry_data: PyNMRSTAR entry object
        :param data_set_id: Data set identifier (bmrb id or filename or user defined id)
        :param auth_tag: Use author sequence numbering True/False default: False
        :return: Chemical shift dictionary {data_set_id:{chain_id:{seq_id:{atom_id:cs_value}},'seq_ids':[1,2,3,4..]}}
        """
        cs_data = {}
        # entry_data = pynmrstar.Entry.from_database(bmrb_id)
        cs_loops = entry_data.get_loops_by_category('Atom_chem_shift')
        logging.debug('Chemical shift loop read')
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
        logging.debug('Getting sequence information from chemical shift loop')
        for cs_id in cs_data.keys():
            for chain in cs_data[cs_id].keys():
                seq_ids = sorted(set([i for i in cs_data[cs_id][chain].keys()]))
                cs_data[cs_id][chain]['seq_ids'] = seq_ids
        return cs_data

    def from_entry_object(self, entry_objects: Union[pynmrstar.Entry, List[pynmrstar.Entry]],
                          auth_tag: bool = False,
                          data_set_id: Union[str, List[str]] = None):
        """
        Extracts chemical shift information one or more PyNMRSTAR entry object

        :param entry_objects: list of pynmrstar entry object  or single entry object
        :param auth_tag: Use author sequence numbering True/False default: False
        :param data_set_id: User defined data set id default: filename
        :return: Chemical shift dictionary {data_set_id:{chain_id:{seq_id:{atom_id:cs_value}},'seq_ids':[1,2,3,4..]}}
        """
        all_cs_data = {}
        if type(entry_objects) is list:
            if type(data_set_id) is list:
                if len(entry_objects) != len(data_set_id):
                    raise ValueError('List of entry objects does not match with data_set_ids')
                else:
                    for ent in entry_objects:
                        data_id = data_set_id[entry_objects.index(ent)]
                        cs_data = self._from_pynmrstar_entry_object(entry_data=ent, data_set_id=data_id,
                                                                    auth_tag=auth_tag)
                        all_cs_data.update(cs_data)
            for ent in entry_objects:
                id = entry_objects.index(ent)
                cs_data = self._from_pynmrstar_entry_object(entry_data=ent, data_set_id='{}_{}'.format(data_set_id, id),
                                                            auth_tag=auth_tag)
                all_cs_data.update(cs_data)
        else:
            if data_set_id is None:
                data_set_id = entry_objects.entry_id
            cs_data = self._from_pynmrstar_entry_object(entry_data=entry_objects, data_set_id=data_set_id,
                                                        auth_tag=auth_tag)
            all_cs_data.update(cs_data)
        return all_cs_data

    def from_file(self, input_file_names: Union[str, List[str]],
                  auth_tag: bool = False,
                  data_set_id: Union[str, List[str]] = None):
        """
        Extracts chemical shift information one or more NMR-STAR files

        :param input_file_names: list of NMR-STAR file names with full path or single file name with full path
        :param auth_tag: Use author sequence numbering True/False default: False
        :param data_set_id: User defined data set id default: filename
        :return: Chemical shift dictionary {data_set_id:{chain_id:{seq_id:{atom_id:cs_value}},'seq_ids':[1,2,3,4..]}}
        """
        if type(input_file_names) is list:
            all_cs_data = {}
            for file_name in input_file_names:
                logging.debug('Reading file {}'.format(file_name))
                if os.path.exists(file_name):
                    entry_data = pynmrstar.Entry.from_file(file_name)
                    data_set_id = os.path.splitext(os.path.basename(file_name))[0]
                    cs_data = self._from_pynmrstar_entry_object(entry_data, data_set_id, auth_tag)
                else:
                    logging.error('File not found {}'.format(file_name))
                    raise IOError('File not found : {}'.format(file_name))
                all_cs_data.update(cs_data)
        else:
            logging.debug('Reading file {}'.format(input_file_names))
            if os.path.exists(input_file_names):
                entry_data = pynmrstar.Entry.from_file(input_file_names)
                if data_set_id is None:
                    data_set_id = os.path.splitext(os.path.basename(input_file_names))[0]
                all_cs_data = self._from_pynmrstar_entry_object(entry_data, data_set_id, auth_tag)
            else:
                logging.error('File not found {}'.format(input_file_names))
                raise IOError('File not found : {}'.format(input_file_names))
        return all_cs_data

    def from_bmrb(self, bmrb_ids: Union[str, List[str]],
                  auth_tag: bool = False):
        """
        Extracts chemical shift information directly from BMRB database for a given BMRB entry or list of entries

        :param bmrb_ids: List of BMRB entries ids or single BMRB ID
        :param auth_tag: Use author sequence numbering True/False default: False
        :return: Chemical shift dictionary {data_set_id:{chain_id:{seq_id:{atom_id:cs_value}},'seq_ids':[1,2,3,4..]}}
        """
        if type(bmrb_ids) is list:
            all_cs_data = {}
            for bmrb_id in bmrb_ids:
                try:
                    logging.debug('Getting entry {} from BMRB'.format(bmrb_id))
                    entry_data = pynmrstar.Entry.from_database(bmrb_id)
                except OSError:
                    entry_data = None
                except KeyError:
                    entry_data = None
                if entry_data is not None:
                    cs_data = self._from_pynmrstar_entry_object(entry_data, bmrb_id, auth_tag)
                else:

                    logging.error('Entry {} not found in public database'.format(bmrb_id))
                    raise IOError('Entry not found in public database: {}'.format(bmrb_id))
                all_cs_data.update(cs_data)
        else:

            try:
                logging.debug('Getting entry {} from BMRB'.format(bmrb_ids))
                entry_data = pynmrstar.Entry.from_database(bmrb_ids)
            except OSError:
                entry_data = None
            except KeyError:
                entry_data = None
            if entry_data is not None:
                all_cs_data = self._from_pynmrstar_entry_object(entry_data, bmrb_ids, auth_tag)
            else:
                logging.error('Entry {} not found in public database'.format(bmrb_ids))
                raise IOError('Entry not found in public database: {}'.format(bmrb_ids))
        return all_cs_data


if __name__ == "__main__":
    p = ChemicalShift().from_bmrb(bmrb_ids='15060')
    print(p)
#     # p.from_file('/Users/kumaran/MyData.str',data_set_id='test')
#     # p.
#

#!/usr/bin/env python
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

_API_URL = "http://api.bmrb.io/v2"
NOTEBOOK = False
_OPACITY = 0.5
_AUTOOPEN = True
__version__ = "2.0.1"

__all__ = ['Spectra', 'Histogram']


# http://webapi.bmrb.wisc.edu/v2/search/chemical_shifts?comp_id=ASP&atom_id=HD2

class Spectra(object):
    """
    Generates HSQC peak positions
    """

    def __init__(self):
        self.oneTOthree = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                           'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                           'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                           'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
        self.threeTOone = {}
        for kk in self.oneTOthree.keys():
            self.threeTOone[self.oneTOthree[kk]] = kk
        if NOTEBOOK:
            plotly.offline.init_notebook_mode(connected=True)

    def get_entry(self, bmrbid=None, filename=None, seq=None, tag='User', nn=3):
        """
        Downloads the chemical shift data for a given entry id or list of entry ids
        :param bmrbid: BMRB ID or list of BMRBIDs
        :param filename: local NMR-STAR filename
        :param seq: amino acids sequence in single letter code
        :param tag: tag associated to seq default 'User'
        :param nn: nearest neighbour effect (3/5/7) default 3
        :return: chemical shift data
        """
        outdata = []
        if bmrbid is None and filename is None and seq is None:
            raise ValueError("BMRB ID or filename or sequence must be given")
        else:
            if seq is not None:
                if nn in [3, 5, 7]:
                    outdata = self.predict_from_seq(seq, tag, nn)
                else:
                    raise ValueError("Wrong nearest neighbour value nn = 3/5/7")
                    outdata = []
            if filename is not None:
                if os.path.exists(filename):
                    indata = pynmrstar.Entry.from_file(filename)
                    fid = os.path.splitext(os.path.basename(filename))[0]
                    cs_data = indata.get_tags(
                        ['_Atom_chem_shift.Comp_index_ID', '_Atom_chem_shift.Comp_ID', '_Atom_chem_shift.Atom_ID',
                         '_Atom_chem_shift.Atom_type', '_Atom_chem_shift.Assigned_chem_shift_list_ID',
                         '_Atom_chem_shift.Val','_Atom_chem_shift.Entity_assembly_ID'])
                    eids = [fid] * len(cs_data['_Atom_chem_shift.Comp_index_ID'])
                    # eids = [fid for i in range(len(cs_data['_Atom_chem_shift.Comp_index_ID']))]
                    eid_cs_data = [eids, cs_data['_Atom_chem_shift.Comp_index_ID'],
                                   cs_data['_Atom_chem_shift.Comp_ID'],
                                   cs_data['_Atom_chem_shift.Atom_ID'],
                                   cs_data['_Atom_chem_shift.Atom_type'],
                                   cs_data['_Atom_chem_shift.Assigned_chem_shift_list_ID'],
                                   cs_data['_Atom_chem_shift.Val'],
                                   cs_data['_Atom_chem_shift.Entity_assembly_ID']]
                    if len(outdata):
                        for i in range(len(eid_cs_data)):
                            outdata[i] = outdata[i] + eid_cs_data[i]
                    else:
                        outdata = eid_cs_data
                else:
                    raise IOError('File not found : {}'.format(filename))
            if bmrbid is not None:
                if type(bmrbid) is list:
                    for eid in bmrbid:
                        try:
                            indata = pynmrstar.Entry.from_database(eid)
                            cs_data = indata.get_tags(
                                ['_Atom_chem_shift.Comp_index_ID', '_Atom_chem_shift.Comp_ID',
                                 '_Atom_chem_shift.Atom_ID',
                                 '_Atom_chem_shift.Atom_type', '_Atom_chem_shift.Assigned_chem_shift_list_ID',
                                 '_Atom_chem_shift.Val','_Atom_chem_shift.Entity_assembly_ID'])
                            eids = [eid] * len(cs_data['_Atom_chem_shift.Comp_index_ID'])
                            # eids = [eid for i in range(len(cs_data['_Atom_chem_shift.Comp_index_ID']))]
                            eid_cs_data = [eids, cs_data['_Atom_chem_shift.Comp_index_ID'],
                                           cs_data['_Atom_chem_shift.Comp_ID'],
                                           cs_data['_Atom_chem_shift.Atom_ID'],
                                           cs_data['_Atom_chem_shift.Atom_type'],
                                           cs_data['_Atom_chem_shift.Assigned_chem_shift_list_ID'],
                                           cs_data['_Atom_chem_shift.Val'],
                                           cs_data['_Atom_chem_shift.Entity_assembly_ID']]
                        except (OSError, IOError) as e:
                            logging.exception("Failed to load entry %s from database with exception %s.", eid, e)
                        if len(outdata):
                            for i in range(len(eid_cs_data)):
                                outdata[i] = outdata[i] + eid_cs_data[i]
                        else:
                            outdata = eid_cs_data
                else:
                    try:
                        indata = pynmrstar.Entry.from_database(bmrbid)
                        cs_data = indata.get_tags(
                            ['_Atom_chem_shift.Comp_index_ID', '_Atom_chem_shift.Comp_ID', '_Atom_chem_shift.Atom_ID',
                             '_Atom_chem_shift.Atom_type', '_Atom_chem_shift.Assigned_chem_shift_list_ID',
                             '_Atom_chem_shift.Val','_Atom_chem_shift.Entity_assembly_ID'])
                        eids = [bmrbid] * len(cs_data['_Atom_chem_shift.Comp_index_ID'])
                        # eids = [bmrbid for i in range(len(cs_data['_Atom_chem_shift.Comp_index_ID']))]
                        eid_cs_data = [eids, cs_data['_Atom_chem_shift.Comp_index_ID'],
                                       cs_data['_Atom_chem_shift.Comp_ID'],
                                       cs_data['_Atom_chem_shift.Atom_ID'],
                                       cs_data['_Atom_chem_shift.Atom_type'],
                                       cs_data['_Atom_chem_shift.Assigned_chem_shift_list_ID'],
                                       cs_data['_Atom_chem_shift.Val'],
                                       cs_data['_Atom_chem_shift.Entity_assembly_ID']]
                        if len(outdata):
                            for i in range(len(eid_cs_data)):
                                outdata[i] = outdata[i] + eid_cs_data[i]
                        else:
                            outdata = eid_cs_data
                    except (OSError, IOError) as e:
                        logging.exception("Failed to load entry %s from database with exception %s.", bmrbid, e)
        return outdata


    @staticmethod
    def _load_pp_dict(atom, nn, filtered=True):
        """
        Internal method to load nearest neighbor statistis from lib folder
        :param atom: atom name
        :param nn: Tri or Pendta or Hepta (3/5/7)
        :param filtered: exclude the outliers
        :return: python dictionary object
        """
        fname = scriptPath + '/data/nn_pp_{}_{}_filtered.txt'.format(atom, nn)
        with open(fname, 'r') as f:
            dat = f.read().split("\n")[:-1]
        if filtered:
            ix = 4
        else:
            ix = 1
        pp_dic = {}
        for l in dat:
            d = l.split("\t")
            pp_dic[d[0]] = float(d[ix])
        return pp_dic

    def predict_from_seq(self, seq=None, tag='User', nn=3):
        """
        Predicts pseudo random coil chemical shifts for a given sequence with standard amino acids
        :param seq: single letter sequence as a string
        :param tag: Name for the sequence
        :param nn: Nearest neighbor effect (3 for Tri-peptide, 5 for Penta and 7 for hepta)
        :return: chemical shift values as a python list
        """
        if seq is not None and nn in [3, 5, 7]:
            atom_list = {'N': 'N', 'H': 'H'}
            nearest_nei = nn
            standards_only = True
            for i in seq:
                if i not in self.oneTOthree.keys():
                    standards_only = False
            if standards_only:
                if nearest_nei == 3:
                    dict_n3 = self._load_pp_dict('N', 1)
                    dict_h3 = self._load_pp_dict('H', 1)
                    dict_n = self._load_pp_dict('N', 0)
                    dict_h = self._load_pp_dict('H', 0)
                    eid = []
                    comp_index_id = []
                    comp_id = []
                    atom_id = []
                    atom_type = []
                    cs_list_id = []
                    val = []
                    for i in [0, len(seq) - 1]:
                        if self.oneTOthree[seq[i]] != "PRO":
                            tp = self.oneTOthree[seq[i]]
                            for atm in atom_list.keys():
                                eid.append(tag)
                                comp_index_id.append(i + 1)
                                comp_id.append(self.oneTOthree[seq[i]])
                                atom_id.append(atm)
                                atom_type.append(atom_list[atm])
                                cs_list_id.append('1')
                                if atm == 'N':
                                    val.append(dict_n[tp])
                                elif atm == 'H':
                                    val.append(dict_h[tp])
                                else:
                                    logging.warning("Something wrong")
                    for i in range(1, len(seq) - 1):
                        if self.oneTOthree[seq[i]] != "PRO":
                            tp3 = '{}-{}-{}'.format(self.oneTOthree[seq[i - 1]], self.oneTOthree[seq[i]],
                                                    self.oneTOthree[seq[i + 1]])
                            for atm in atom_list.keys():
                                eid.append(tag)
                                comp_index_id.append(i + 1)
                                comp_id.append(self.oneTOthree[seq[i]])
                                atom_id.append(atm)
                                atom_type.append(atom_list[atm])
                                cs_list_id.append('1')
                                if atm == 'N':
                                    val.append(dict_n3[tp3])
                                elif atm == 'H':
                                    val.append(dict_h3[tp3])
                                else:
                                    logging.warning("Something wrong")
                elif nearest_nei == 5:
                    dict_n5 = self._load_pp_dict('N', 2)
                    dict_h5 = self._load_pp_dict('H', 2)
                    dict_n3 = self._load_pp_dict('N', 1)
                    dict_h3 = self._load_pp_dict('H', 1)
                    dict_n = self._load_pp_dict('N', 0)
                    dict_h = self._load_pp_dict('H', 0)
                    eid = []
                    comp_index_id = []
                    comp_id = []
                    atom_id = []
                    atom_type = []
                    cs_list_id = []
                    val = []
                    for i in [0, 1, len(seq) - 2, len(seq) - 1]:
                        if self.oneTOthree[seq[i]] != "PRO":
                            if i == 0 or i == len(seq) - 1:
                                tp = self.oneTOthree[seq[i]]
                                for atm in atom_list.keys():
                                    eid.append(tag)
                                    comp_index_id.append(i + 1)
                                    comp_id.append(self.oneTOthree[seq[i]])
                                    atom_id.append(atm)
                                    atom_type.append(atom_list[atm])
                                    cs_list_id.append('1')
                                    if atm == 'N':
                                        val.append(dict_n[tp])
                                    elif atm == 'H':
                                        val.append(dict_h[tp])
                                    else:
                                        logging.warning("Something wrong")
                            else:
                                tp3 = '{}-{}-{}'.format(self.oneTOthree[seq[i - 1]], self.oneTOthree[seq[i]],
                                                        self.oneTOthree[seq[i + 1]])
                                for atm in atom_list.keys():
                                    eid.append(tag)
                                    comp_index_id.append(i + 1)
                                    comp_id.append(self.oneTOthree[seq[i]])
                                    atom_id.append(atm)
                                    atom_type.append(atom_list[atm])
                                    cs_list_id.append('1')
                                    if atm == 'N':
                                        val.append(dict_n3[tp3])
                                    elif atm == 'H':
                                        val.append(dict_h3[tp3])
                                    else:
                                        logging.warning("Something wrong")
                    for i in range(2, len(seq) - 2):
                        if self.oneTOthree[seq[i]] != "PRO":
                            tp5 = '{}-{}-{}-{}-{}'.format(self.oneTOthree[seq[i - 2]], self.oneTOthree[seq[i - 1]],
                                                          self.oneTOthree[seq[i]],
                                                          self.oneTOthree[seq[i + 1]], self.oneTOthree[seq[i + 2]])
                            tp3 = '{}-{}-{}'.format(self.oneTOthree[seq[i - 1]], self.oneTOthree[seq[i]],
                                                    self.oneTOthree[seq[i + 1]])
                            for atm in atom_list.keys():
                                eid.append(tag)
                                comp_index_id.append(i + 1)
                                comp_id.append(self.oneTOthree[seq[i]])
                                atom_id.append(atm)
                                atom_type.append(atom_list[atm])
                                cs_list_id.append('1')
                                if atm == 'N':
                                    try:
                                        val.append(dict_n5[tp5])
                                    except KeyError:
                                        val.append(dict_n3[tp3])
                                elif atm == 'H':
                                    try:
                                        val.append(dict_h5[tp5])
                                    except KeyError:
                                        val.append(dict_h3[tp3])
                                else:
                                    logging.warning("Something wrong")
                else:
                    dict_n7 = self._load_pp_dict('N', 3)
                    dict_h7 = self._load_pp_dict('H', 3)
                    dict_n5 = self._load_pp_dict('N', 2)
                    dict_h5 = self._load_pp_dict('H', 2)
                    dict_n3 = self._load_pp_dict('N', 1)
                    dict_h3 = self._load_pp_dict('H', 1)
                    dict_n = self._load_pp_dict('N', 0)
                    dict_h = self._load_pp_dict('H', 0)
                    eid = []
                    comp_index_id = []
                    comp_id = []
                    atom_id = []
                    atom_type = []
                    cs_list_id = []
                    val = []
                    for i in [0, 1, 2, len(seq) - 3, len(seq) - 2, len(seq) - 1]:
                        if self.oneTOthree[seq[i]] != "PRO":
                            if i == 0 or i == len(seq) - 1:
                                tp = self.oneTOthree[seq[i]]
                                for atm in atom_list.keys():
                                    eid.append(tag)
                                    comp_index_id.append(i + 1)
                                    comp_id.append(self.oneTOthree[seq[i]])
                                    atom_id.append(atm)
                                    atom_type.append(atom_list[atm])
                                    cs_list_id.append('1')
                                    if atm == 'N':
                                        val.append(dict_n[tp])
                                    elif atm == 'H':
                                        val.append(dict_h[tp])
                                    else:
                                        logging.warning("Something wrong")
                            elif i == 1 or i == len(seq) - 2:
                                tp3 = '{}-{}-{}'.format(self.oneTOthree[seq[i - 1]], self.oneTOthree[seq[i]],
                                                        self.oneTOthree[seq[i + 1]])
                                for atm in atom_list.keys():
                                    eid.append(tag)
                                    comp_index_id.append(i + 1)
                                    comp_id.append(self.oneTOthree[seq[i]])
                                    atom_id.append(atm)
                                    atom_type.append(atom_list[atm])
                                    cs_list_id.append('1')
                                    if atm == 'N':
                                        val.append(dict_n3[tp3])
                                    elif atm == 'H':
                                        val.append(dict_h3[tp3])
                                    else:
                                        logging.warning("Something wrong")
                            else:
                                tp5 = '{}-{}-{}-{}-{}'.format(self.oneTOthree[seq[i - 2]], self.oneTOthree[seq[i - 1]],
                                                              self.oneTOthree[seq[i]],
                                                              self.oneTOthree[seq[i + 1]], self.oneTOthree[seq[i + 2]])
                                tp3 = '{}-{}-{}'.format(self.oneTOthree[seq[i - 1]], self.oneTOthree[seq[i]],
                                                        self.oneTOthree[seq[i + 1]])
                                for atm in atom_list.keys():
                                    eid.append(tag)
                                    comp_index_id.append(i + 1)
                                    comp_id.append(self.oneTOthree[seq[i]])
                                    atom_id.append(atm)
                                    atom_type.append(atom_list[atm])
                                    cs_list_id.append('1')
                                    if atm == 'N':
                                        try:
                                            val.append(dict_n5[tp5])
                                        except KeyError:
                                            val.append(dict_n3[tp3])
                                    elif atm == 'H':
                                        try:
                                            val.append(dict_h5[tp5])
                                        except KeyError:
                                            val.append(dict_h3[tp3])
                                    else:
                                        logging.warning("Something wrong")
                    for i in range(3, len(seq) - 3):
                        if self.oneTOthree[seq[i]] != "PRO":
                            tp7 = '{}-{}-{}-{}-{}-{}-{}'.format(self.oneTOthree[seq[i - 3]],
                                                                self.oneTOthree[seq[i - 2]],
                                                                self.oneTOthree[seq[i - 1]], self.oneTOthree[seq[i]],
                                                                self.oneTOthree[seq[i + 1]],
                                                                self.oneTOthree[seq[i + 2]],
                                                                self.oneTOthree[seq[i + 3]])
                            tp5 = '{}-{}-{}-{}-{}'.format(self.oneTOthree[seq[i - 2]], self.oneTOthree[seq[i - 1]],
                                                          self.oneTOthree[seq[i]],
                                                          self.oneTOthree[seq[i + 1]], self.oneTOthree[seq[i + 2]])
                            tp3 = '{}-{}-{}'.format(self.oneTOthree[seq[i - 1]], self.oneTOthree[seq[i]],
                                                    self.oneTOthree[seq[i + 1]])
                            for atm in atom_list.keys():
                                eid.append(tag)
                                comp_index_id.append(i + 1)
                                comp_id.append(self.oneTOthree[seq[i]])
                                atom_id.append(atm)
                                atom_type.append(atom_list[atm])
                                cs_list_id.append('1')
                                if atm == 'N':
                                    try:
                                        val.append(dict_n7[tp7])
                                    except KeyError:
                                        try:
                                            val.append(dict_n5[tp5])
                                        except KeyError:
                                            val.append(dict_n3[tp3])
                                elif atm == 'H':
                                    try:
                                        val.append(dict_h7[tp7])
                                    except KeyError:
                                        try:
                                            val.append(dict_h5[tp5])
                                        except KeyError:
                                            val.append(dict_h3[tp3])
                                else:
                                    logging.warning("Something wrong")

                out_data = [eid, comp_index_id, comp_id, atom_id, atom_type, cs_list_id, val]
            else:
                logging.exception("Sequence contains nonstandard amino acid code. Random coil prediction is not possible with non "
                      "standard amino acids")
                out_data = []
        else:
            out_data = []
        return out_data


    @staticmethod
    def check_hsqc_peaks(hsqcdata):
        hsqcnew = [[],[],[],[],[],[]]
        logging.debug(hsqcdata)
        for i in range(len(hsqcdata[0])):
            try:
                hcs = float(hsqcdata[1][i])
                ncs = float(hsqcdata[2][i])
                hsqcnew[0].append(hsqcdata[0][i])
                hsqcnew[1].append(hsqcdata[1][i])
                hsqcnew[2].append(hsqcdata[2][i])
                hsqcnew[3].append(hsqcdata[3][i])
                hsqcnew[4].append(hsqcdata[4][i])
            except TypeError:
                pass
        return hsqcnew


    @staticmethod
    def convert_to_n15hsqc_peaks(csdata):
        """
        Converts the output from get_entry into hsqc peak positions
        :param csdata: output from get_entry
        :return: easy to plot hsqc peak positions
        """

        sidechainres = ['ARG', 'GLN', 'ASN', 'HIS', 'TRP', 'LYS']
        sidechains = {
            'ARG-HH11': ['HH11''NH1'],
            'ARG-HH12': ['HH12', 'NH1'],
            'ARG-HH21': ['HH21', 'NH2'],
            'ARG-HH22': ['HH22', 'NH2'],
            'ARG-HE': ['HE', 'NE'],
            'GLN-HE21': ['HE21', 'NE2'],
            'GLN-HE22': ['HE22', 'NE2'],
            'ASN-HD21': ['HD21', 'ND2'],
            'ASN-HD22': ['HD22', 'ND2'],
            'HIS-HD1': ['HD1''ND1'],
            'HIS-HE2': ['HE2', 'NE2'],
            'TRP-HE1': ['HE1', 'NE1'],
            'LYS-HZ': ['HZ', 'NZ'],
            'LYS-HZ1': ['HZ1', 'NZ'],
            'LYS-HZ2': ['HZ2', 'NZ'],
            'LYS-HZ3': ['HZ3', 'NZ']
        }
        outdata = [[], [], [], [], []]
        for i in range(len(csdata[0])):
            atomid = '{}-{}-{}-{}-{}'.format(csdata[0][i], csdata[1][i],  csdata[2][i], csdata[5][i], csdata[7][i])
            if csdata[3][i] == "H":
                if atomid not in outdata[0]:
                    outdata[0].append(atomid)
                    outdata[1].append(csdata[6][i])
                    outdata[2].append(None)
                    outdata[3].append(csdata[3][i])
                    outdata[4].append(None)
                else:
                    outdata[1][outdata[0].index(atomid)] = csdata[6][i]
                    outdata[3][outdata[0].index(atomid)] = csdata[3][i]
            if csdata[3][i] == "N":
                if atomid not in outdata[0]:
                    outdata[0].append(atomid)
                    outdata[2].append(csdata[6][i])
                    outdata[1].append(None)
                    outdata[4].append(csdata[3][i])
                    outdata[3].append(None)
                else:
                    outdata[2][outdata[0].index(atomid)] = csdata[6][i]
                    outdata[4][outdata[0].index(atomid)] = (csdata[3][i])
            if csdata[2][i] in sidechainres:
                for k in sidechains.keys():
                    if k.split("-")[0] == csdata[2][i]:
                        atomid = '{}-{}-{}-{}-{}'.format(csdata[0][i], csdata[1][i], k, csdata[5][i],csdata[7][i])
                        if csdata[3][i] in sidechains[k] and csdata[4][i] == "H":
                            if atomid not in outdata[0]:
                                outdata[0].append(atomid)
                                outdata[1].append(csdata[6][i])
                                outdata[2].append(None)
                                outdata[3].append(csdata[3][i])
                                outdata[4].append(None)

                            else:
                                outdata[1][outdata[0].index(atomid)] = csdata[6][i]
                                outdata[3][outdata[0].index(atomid)] = csdata[3][i]
                        if csdata[3][i] in sidechains[k] and csdata[4][i] == "N":
                            if atomid not in outdata[0]:
                                outdata[0].append(atomid)
                                outdata[2].append(csdata[6][i])
                                outdata[1].append(None)
                                outdata[4].append(csdata[3][i])
                                outdata[3].append(None)
                            else:
                                outdata[2][outdata[0].index(atomid)] = csdata[6][i]
                                outdata[4][outdata[0].index(atomid)] = (csdata[3][i])
        return outdata


    def n15hsqc(self, bmrbid=None, filename=None, seq=None, nn=3, colorby=None, groupbyres=False,
                outfilename='n15hsqc.html'):
        """
        Plots hsqc peak positions for a given list of BMRB ids
        :param bmrbid: entry id or list of entry ids
        :param seq: sequence in one letter code
        :param filename: local NMR-STAR filename
        :param nn: 3/5/7 for Tri or Penta or Hepta peptide model prediction
        :param colorby: Color by res/entry (default res for single entry/ entry for multiple entries)
        :param groupbyres: if TRUE connects the same seq ids by line; default False
        :param outfilename: Output filename
        :return: plotly plot object or html file
        """
        if nn == 3:
            tag = "TriPeptide"
        elif nn == 5:
            tag = "PentaPeptide"
        elif nn == 7:
            tag = "HeptaPeptide"
        else:
            tag = "Not a valid model"
        title = 'Simulated <sup>1</sup>H-<sup>15</sup>N  HSQC peak positions'

        csdata = self.get_entry(bmrbid, filename, seq, tag, nn)

        if len(csdata):
            hsqcdata = self.check_hsqc_peaks(self.convert_to_n15hsqc_peaks(csdata))
        else:
            hsqcdata = []

        if colorby == 'entry':
            idx = 0
        elif colorby == 'res':
            idx = 2
        elif type(bmrbid) is list and len(bmrbid) > 1:
            idx = 0
        elif sum(1 for i in [bmrbid, filename, seq] if i is not None) > 1:
            idx = 0
        else:
            idx = 2

        if len(hsqcdata):
            groups = set([k.split("-")[idx] for k in hsqcdata[0]])
        else:
            groups = []
        data_sets = {}
        for gid in groups:
            data_sets[gid] = [[], [], [], []]
            for i in range(len(hsqcdata[0])):
                if hsqcdata[0][i].split("-")[idx] == gid:
                    data_sets[gid][0].append(hsqcdata[1][i])
                    data_sets[gid][1].append(hsqcdata[2][i])
                    data_sets[gid][2].append(hsqcdata[0][i])
                    try:
                        data_sets[gid][3].append(
                            '{}{}'.format(hsqcdata[0][i].split("-")[1], self.threeTOone[hsqcdata[0][i].split("-")[2]]))
                    except KeyError:
                        data_sets[gid][3].append(
                            '{}{}'.format(hsqcdata[0][i].split("-")[1], hsqcdata[0][i].split("-")[2]))

        if groupbyres:
            groups2 = set(["-".join(k.split("-")[1:4]) for k in hsqcdata[0]])
            data_sets2 = {}
            for gid in groups2:
                data_sets2[gid] = [[], [], [], []]
                for i in range(len(hsqcdata[0])):
                    if "-".join(hsqcdata[0][i].split("-")[1:4]) == gid:
                        data_sets2[gid][0].append(hsqcdata[1][i])
                        data_sets2[gid][1].append(hsqcdata[2][i])
                        data_sets2[gid][2].append(hsqcdata[0][i])
                        try:
                            one_letter_code = self.threeTOone[hsqcdata[0][i].split("-")[2]]
                        except KeyError:
                            one_letter_code = hsqcdata[0][i].split("-")[2]
                        seq_id = hsqcdata[0][i].split("-")[1]
                        data_sets2[gid][3].append('{}{}'.format(seq_id, one_letter_code))

        data = []
        for k in data_sets.keys():
            data.append(plotly.graph_objs.Scatter(x=data_sets[k][0],
                                                  y=data_sets[k][1],
                                                  text=data_sets[k][2],
                                                  mode='markers',
                                                  opacity=_OPACITY,
                                                  name=k))

        if groupbyres:
            for k in data_sets2.keys():
                data.append(plotly.graph_objs.Scatter(x=data_sets2[k][0],
                                                      y=data_sets2[k][1],
                                                      text=data_sets2[k][2],
                                                      mode='lines',
                                                      name=k,
                                                      opacity=_OPACITY,
                                                      showlegend=False)
                            )

        layout = plotly.graph_objs.Layout(
            xaxis=dict(autorange='reversed',
                       title='<sup>1</sup>H (ppm)',
                       type='linear'),
            yaxis=dict(autorange='reversed',
                       title='<sup>15</sup>N (ppm)',
                       type='linear'),
            showlegend=True,
            hovermode='closest',
            title=title)
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        if len(data):
            if NOTEBOOK:
                plotly.offline.iplot(fig)
            else:
                plotly.offline.plot(fig, filename=outfilename, auto_open=_AUTOOPEN)
                logging.info("Output written on {}".format(outfilename))
            return True
        else:
            raise ValueError("No amide proton nitrogen chemical shifts found")
            return False



class Histogram(object):
    """
    Generates checmial shift histograms from BMRB database
    """

    def __init__(self):
        if NOTEBOOK:
            plotly.offline.init_notebook_mode(connected=True)

    @staticmethod
    def get_histogram_api(residue, atom, filtered=True, sd_limit=10, normalized=False, ambiguity = '*'):
        """
        Downloads chemical shift data for a given atom in a residue using BMRB API
        :param residue: Residue name in standard 3 letter code
        :param atom: IUIPAC atom name
        :param filtered: True/False Filters based on standard deviation cutoff Default:True
        :param sd_limit: Number of time Standard deviation for filtering default: 10
        :param normalized: True/False Plots either Count/Density default: False
        :param ambiguity: NMR-STAR Ambiguity code default:'*'
        :return: Plotly object
        """
        standard = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS',
                    'ASP', 'SER', 'LYS', 'PRO', 'ASN',
                    'VAL', 'THR', 'HIS', 'TRP', 'PHE',
                    'ALA', 'MET', 'LEU', 'ARG', 'TYR',
                    'A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']

        if residue == "*" and atom == "*":
            raise ValueError("Getting full database will overload the memory! Please chose a residue or atom.")
        elif residue == "*":
            url = Request(_API_URL + "/search/chemical_shifts?atom_id={}".format(atom))
        elif atom == "*":
            url = Request(_API_URL + "/search/chemical_shifts?comp_id={}".format(residue))
        else:
            url = Request(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom))
        url.add_header('Application', 'PyBMRB')
        r = urlopen(url)
        dump = json.loads(r.read())
        if residue == "*":
            d = [i for i in dump['data'] if i[dump['columns'].index('Atom_chem_shift.Comp_ID')] in standard]
        else:
            d = [i for i in dump['data']]

        alist = set(['{}-{}'.format(i[dump['columns'].index('Atom_chem_shift.Comp_ID')],
                                    i[dump['columns'].index('Atom_chem_shift.Atom_ID')]) for i in d])
        if len(alist) > 1:
            data = []
            for atm in alist:
                res = atm.split("-")[0]
                at = atm.split("-")[1]
                if ambiguity == '*':
                    x = [i[dump['columns'].index('Atom_chem_shift.Val')] for i in d
                         if i[dump['columns'].index('Atom_chem_shift.Comp_ID')] == res and
                         i[dump['columns'].index('Atom_chem_shift.Atom_ID')] == at]
                else:
                    x = [i[dump['columns'].index('Atom_chem_shift.Val')] for i in d
                         if i[dump['columns'].index('Atom_chem_shift.Comp_ID')] == res and
                         i[dump['columns'].index('Atom_chem_shift.Atom_ID')] == at and
                         i[dump['columns'].index('Atom_chem_shift.Ambiguity_code')] == ambiguity]

                if filtered:
                    mean = np.mean(x)
                    sd = np.std(x)
                    lb = mean - (sd_limit * sd)
                    ub = mean + (sd_limit * sd)
                    x = [i for i in x if lb < i < ub]
                if len(x) == 0:
                    logging.warning('{} has no data at BMRB. Please check the atom nomenclature.'.format(atm))
                if normalized:
                    data.append(plotly.graph_objs.Histogram(x=x, name=atm,
                                                            histnorm='probability', opacity=_OPACITY))

                else:
                    data.append(plotly.graph_objs.Histogram(x=x, name=atm, opacity=_OPACITY))

        else:
            if ambiguity == '*':
                x = [i[dump['columns'].index('Atom_chem_shift.Val')] for i in d]
            else:
                x = [i[dump['columns'].index('Atom_chem_shift.Val')] for i in d
                     if i[dump['columns'].index('Atom_chem_shift.Ambiguity_code')] == ambiguity]
            if len(x) == 0:
                logging.warning('{}-{} has no data at BMRB. Please check the atom nomenclature.'.format(residue, atom))
            if filtered:
                mean = np.mean(x)
                sd = np.std(x)
                lb = mean - (sd_limit * sd)
                ub = mean + (sd_limit * sd)
                x = [i for i in x if lb < i < ub]

            if normalized:

                data = [plotly.graph_objs.Histogram(x=x, name="{}-{}".format(residue, atom),
                                                    histnorm='probability', opacity=_OPACITY)]
            else:
                data = [plotly.graph_objs.Histogram(x=x, name="{}-{}".format(residue, atom), opacity=_OPACITY)]
        return data

    @staticmethod
    def get_conditional_histogram_api(residue, atom, atomlist, cslist, filtered=True,
                                      sd_limit=10, normalized=False, ambiguity='*'):
        """
        Downloads the chemical shift data for a given atom in a residue and filters based on a other atom chemical
        shifts in the same residue
        :param residue: Residue name in standard 3 letter code
        :param atom: IUIPAC atom name
        :param atomlist: atom list in IUPAC format as a list
        :param cslist: corresponding chemical shift list
        :param filtered: True/False Filters based on standard deviation cutoff Default:True
        :param sd_limit:  Number of time Standard deviation for filtering default: 10
        :param normalized: True/False Plots either Count/Density default: False
        :param ambiguity: NMR-STAR Ambiguity code default:'*'
        :return: Plotly object
        """
        url = Request(_API_URL + "/search/chemical_shifts?comp_id={}".format(residue))
        url.add_header('Application', 'PyBMRB')
        r = urlopen(url)
        dump = json.loads(r.read())
        if ambiguity == '*':
            d1=dump
        else:
            d1={}
            d1['columns']=dump['columns']
            d1['data'] = [i for i in dump['data'] if i[dump['columns'].index('Atom_chem_shift.Ambiguity_code')] == ambiguity]
        d = {}
        entry_id_index = d1['columns'].index('Atom_chem_shift.Entry_ID')
        seq_id_index = d1['columns'].index('Atom_chem_shift.Comp_index_ID')
        res_id_index = d1['columns'].index('Atom_chem_shift.Comp_ID')
        atom_id_index = d1['columns'].index('Atom_chem_shift.Atom_ID')
        cs_id_index = d1['columns'].index('Atom_chem_shift.Val')
        for i in d1['data']:
            entry_id = i[entry_id_index]
            seq_id = i[seq_id_index]
            res_id = i[res_id_index]
            atom_id = i[atom_id_index]
            d['{}-{}-{}-{}'.format(entry_id, seq_id, res_id, atom_id)] = i[cs_id_index]
        filter_list = []
        for k in d.keys():
            for i in range(len(atomlist)):
                if 'H' in atomlist[i]:
                    epsilon = 0.1
                else:
                    epsilon = 0.5
                if k.split("-")[-1] == atomlist[i] and (cslist[i] + epsilon < d[k] or d[k] < cslist[i] - epsilon):
                    filter_list.append('{}-{}'.format(k.split("-")[0], k.split("-")[1]))
        x = []
        filter_list = list(set(filter_list))
        for k in d.keys():
            atm_id = '{}-{}'.format(k.split("-")[0], k.split("-")[1])
            if atm_id not in filter_list and k.split("-")[2] == residue and k.split("-")[3] == atom:
                x.append(d[k])
        if filtered:
            mean = np.mean(x)
            sd = np.std(x)
            lb = mean - (sd_limit * sd)
            ub = mean + (sd_limit * sd)
            x = [i for i in x if lb < i < ub]
        filter_values = ''
        for i in range(len(atomlist)):
            if i == len(atomlist) - 1:
                filter_values += '{}:{}'.format(atomlist[i], cslist[i])
            else:
                filter_values += '{}:{},'.format(atomlist[i], cslist[i])
        if normalized:
            data = plotly.graph_objs.Histogram(x=x, name="{}-{}({})".format(residue, atom, filter_values),
                                               histnorm='probability', opacity=_OPACITY)
        else:
            data = plotly.graph_objs.Histogram(x=x,
                                               name="{}-{}({})".format(residue, atom, filter_values, opacity=_OPACITY))
        return data

    @staticmethod
    def get_histogram2d_api(residue, atom1, atom2, filtered=True, sd_limit=10, normalized=False,
                            ambiguity1='*', ambiguity2='*'):
        """
        Calculates the correlation between two atoms in the same residue
        :param residue: Residue name in standard 3 letter code
        :param atom1: IUIPAC atom name
        :param atom2: IUIPAC atom name
        :param filtered:  True/False Filters based on standard deviation cutoff Default:True
        :param sd_limit: Number of time Standard deviation for filtering default: 10
        :param normalized: True/False Plots either Count/Density default: False
        :param ambiguity1: NMR-STAR Ambiguity for atom1 code default:'*'
        :param ambiguity2: NMR-STAR Ambiguity for atom2 code default:'*'
        :return: Plotly object
        """
        url1 = Request(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom1))
        url2 = Request(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom2))
        url1.add_header('Application', 'PyBMRB')
        url2.add_header('Application', 'PyBMRB')
        r1 = urlopen(url1)
        r2 = urlopen(url2)
        dump1 = json.loads(r1.read())
        if ambiguity1 == '*':
            d1=dump1
        else:
            d1={}
            d1['columns']=dump1['columns']
            d1['data'] = [i for i in dump1['data'] if i[dump1['columns'].index('Atom_chem_shift.Ambiguity_code')] == ambiguity1]
        d = {}
        for i in d1['data']:
            entry_id = i[d1['columns'].index('Atom_chem_shift.Entry_ID')]
            seq_id = i[d1['columns'].index('Atom_chem_shift.Comp_index_ID')]
            d["{}-{}".format(entry_id, seq_id)] = i[d1['columns'].index('Atom_chem_shift.Val')]
        # x = [i[d1['columns'].index('Atom_chem_shift.Val')] for i in d1['data']]
        dump2 = json.loads(r2.read())
        if ambiguity2 == '*':
            d2=dump2
        else:
            d2={}
            d2['columns'] = dump2['columns']
            d2['data'] = [i for i in dump2['data'] if i[dump2['columns'].index('Atom_chem_shift.Ambiguity_code')] == ambiguity2]
        x = []
        y = []
        for i in d2['data']:
            entry_id = i[d2['columns'].index('Atom_chem_shift.Entry_ID')]
            seq_id = i[d2['columns'].index('Atom_chem_shift.Comp_index_ID')]
            try:
                k = "{}-{}".format(entry_id, seq_id)
                x.append(d[k])
                y.append(i[d2['columns'].index('Atom_chem_shift.Val')])
            except KeyError:
                pass

        if len(x) == 0:
            raise ValueError("No data found for {}".format(atom2))
        if len(y) == 0:
            raise ValueError("No data found for {}".format(atom1))

        # y = [i[d2['columns'].index('Atom_chem_shift.Val')] for i in d2['data']]
        if filtered:
            meanx = np.mean(x)
            meany = np.mean(y)
            sdx = np.std(x)
            sdy = np.std(y)
            lbx = meanx - (sd_limit * sdx)
            lby = meany - (sd_limit * sdy)
            ubx = meanx + (sd_limit * sdx)
            uby = meany + (sd_limit * sdy)
            x1 = [x[i] for i in range(len(x)) if lbx < x[i] < ubx and lby < y[i] < uby]
            y1 = [y[i] for i in range(len(y)) if lby < y[i] < uby and lbx < x[i] < ubx]
            x = x1
            y = y1
        if 'H' in atom2:
            binsizey = 0.05
        else:
            binsizey = 0.25

        if 'H' in atom1:
            binsizex = 0.05
        else:
            binsizex = 0.25
        nbinsx = int(round((max(x) - min(x)) / binsizex))
        nbinsy = int(round((max(y) - min(y)) / binsizey))
        if normalized:
            data = [plotly.graph_objs.Histogram2dContour(x=x, y=y, nbinsy=nbinsy, nbinsx=nbinsx,
                                                         histnorm='probability', colorscale='Jet'),
                    plotly.graph_objs.Histogram(
                        y=y,
                        xaxis='x2',
                        name="{}-{}".format(residue, atom2),
                        nbinsy=nbinsx,
                        histnorm='probability'
                    ),
                    plotly.graph_objs.Histogram(
                        x=x,
                        yaxis='y2',
                        nbinsx=nbinsy,
                        name="{}-{}".format(residue, atom1),
                        histnorm='probability'
                    )
                    ]
        else:
            data = [plotly.graph_objs.Histogram2dContour(x=x, y=y, colorscale='Jet'),
                    plotly.graph_objs.Histogram(
                        y=y,
                        xaxis='x2',
                        nbinsy=nbinsy,
                        name="{}-{}".format(residue, atom2)
                    ),
                    plotly.graph_objs.Histogram(
                        x=x,
                        nbinsx=nbinsx,
                        yaxis='y2',
                        name="{}-{}".format(residue, atom1)
                    )
                    ]
        return data

    def hist(self, residue=None, atom=None, atom_list=None, filtered=True, sd_limit=10,
             normalized=False, ambiguity = '*', outfilename='hist.html'):
        """
            Chemical shift histogram from BMRB database
        :param residue: 3 letter amino acid code
        :param atom: IUPAC atom name
        :param atom_list: list of atoms example: ['ALA-CA','GLY-N']
        :param filtered: True/False Filters based on standard deviation cutoff Default:True
        :param sd_limit:  Number of time Standard deviation for filtering default: 10
        :param normalized: True/False Plots either Count/Density default: False
        :param ambiguity: NMR-STAR Ambiguity code default:'*'
        :param outfilename: output file name
        :return: writes output in a html file
        """
        logging.info("This may take time to gather data from entire BMRB database!")
        if normalized:
            count = 'Density'
        else:
            count = 'Count'
        if residue is None and atom is None and atom_list is None:
            raise ValueError('At least one of the three parameters should be given: residue, atom or atomlist')
        elif (residue is not None and atom is not None and atom_list is not None) or \
                (residue is None and atom is not None and atom_list is not None) or \
                (residue is not None and atom is None and atom_list is not None):

            if residue is None:
                residue = "*"
            if atom is None:
                atom = "*"
            d1 = self.get_histogram_api(residue, atom, filtered, sd_limit, normalized,ambiguity)
            d2 = []
            for atm in atom_list:
                r = atm.split("-")[0]
                a = atm.split("-")[1]
                for d in self.get_histogram_api(r, a, filtered, sd_limit, normalized,ambiguity):
                    d2.append(d)
            data = d1 + d2
            layout = plotly.graph_objs.Layout(
                barmode='overlay',
                xaxis=dict(autorange='reversed', title='Chemical Shift [ppm]'),
                yaxis=dict(title=count))
            fig = plotly.graph_objs.Figure(data=data, layout=layout)
            if NOTEBOOK:
                plotly.offline.iplot(fig)
            else:
                plotly.offline.plot(fig, filename=outfilename, auto_open=_AUTOOPEN)
                logging.info("Output written on {}".format(outfilename))
        elif residue is None and atom is None and atom_list is not None:
            self.multiple_atom(atom_list, filtered, sd_limit, normalized,ambiguity,outfilename,)
        elif residue is not None and atom is not None and atom_list is None:
            self.single_atom(residue, atom, filtered, sd_limit, normalized, ambiguity, outfilename)
        elif residue is not None and atom is None and atom_list is None:
            self.single_atom(residue, '*', filtered, sd_limit, normalized, ambiguity, outfilename)
        elif residue is None and atom is not None and atom_list is None:
            self.single_atom('*', atom, filtered, sd_limit, normalized, ambiguity, outfilename)
        else:
            raise ValueError("Not a valid option")

    def hist2d(self, residue, atom1, atom2, filtered=True, sd_limit=10, normalized=False,
               ambiguity1='*', ambiguity2='*',outfilename='hist2d.html'):
        """
        Generates chemical shift correlation plots for a given two atoms in a given amino acid
        :param residue: 3 letter amino acid code
        :param atom1: IUPAC atom name
        :param atom2: IUPAC atom name
        :param filtered:True/False Filters based on standard deviation cutoff Default:True
        :param sd_limit: Number of time Standard deviation for filtering default: 10
        :param normalized: True/False Plots either Count/Density default: False
        :param outfilename: output file name
        :param ambiguity1: NMR-STAR Ambiguity code for atom1 default:'*'
        :param ambiguity2: NMR-STAR Ambiguity code for atom2 default:'*'
        :return: writes output in a html file
        """
        logging.info("This may take time to gather data from entire BMRB database!")
        layout = plotly.graph_objs.Layout(
            autosize=True,
            xaxis=dict(
                autorange='reversed',
                zeroline=False,
                domain=[0, 0.85],
                showgrid=True,
                title='{}-{} [ppm]'.format(residue, atom1)

            ),
            yaxis=dict(
                autorange='reversed',
                zeroline=False,
                domain=[0, 0.85],
                showgrid=True,
                title='{}-{} [ppm]'.format(residue, atom2)

            ),
            xaxis2=dict(
                zeroline=False,
                domain=[0.85, 1],
                showgrid=True

            ),
            yaxis2=dict(
                zeroline=False,
                domain=[0.85, 1],
                showgrid=True
            ),
            hovermode='closest',
            showlegend=False
        )
        data = self.get_histogram2d_api(residue, atom1, atom2, filtered, sd_limit, normalized, ambiguity1, ambiguity2)
        fig = plotly.graph_objs.Figure(data=data, layout=layout)

        if NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=outfilename, auto_open=_AUTOOPEN)
            logging.info("Output written on {}".format(outfilename))

    def single_atom(self, residue, atom, filtered=True, sd_limit=10, normalized=False, ambiguity = '*',outfilename='hist.html'):
        """
        Generates histgram for a given atom in a given amino acid
        :param residue: 3 letter amino acid code
        :param atom: IUPAC atom name
        :param filtered: True/False Filters based on standard deviation cutoff Default:True
        :param sd_limit: Number of time Standard deviation for filtering default: 10
        :param normalized: True/False Plots either Count/Density default: False
        :param outfilename: output file name
        :return: writes output in a html file
        """
        if normalized:
            count = 'Density'
        else:
            count = 'Count'
        layout = plotly.graph_objs.Layout(
            barmode='overlay',
            xaxis=dict(autorange='reversed', title='Chemical Shift [ppm]'),
            yaxis=dict(title=count))
        data = self.get_histogram_api(residue, atom, filtered, sd_limit, normalized, ambiguity)
        fig = plotly.graph_objs.Figure(data=data, layout=layout)

        if NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=outfilename, auto_open=_AUTOOPEN)
            logging.info("Output written on {}".format(outfilename))

    def multiple_atom(self, atom_list, filtered=True, sd_limit=10, normalized=False, ambiguity = '*',outfilename='hist.html'):
        """
        Generates histogram for a given list of atoms from various amino acids
        :param atom_list: atom list example ['ALA:CA','GLY:CA','ALA:HA']
        :param filtered: True/False Filters based on standard deviation cutoff Default:True
        :param sd_limit: Number of time Standard deviation for filtering default: 10
        :param normalized: True/False Plots either Count/Density default: False
        :param ambiguity: NMR-STAR Ambiguity code default:'*'
        :param outfilename: output file name
        :return: writes output in a html file
        """
        if normalized:
            count = 'Density'
        else:
            count = 'Count'
        data = []
        for atm in atom_list:
            residue = atm.split("-")[0]
            atom = atm.split("-")[1]
            for dd in self.get_histogram_api(residue, atom, filtered, sd_limit, normalized,ambiguity):
                data.append(dd)
        layout = plotly.graph_objs.Layout(
            barmode='overlay',
            xaxis=dict(autorange='reversed', title='Chemical Shift [ppm]'),
            yaxis=dict(title=count))
        fig = plotly.graph_objs.Figure(data=data, layout=layout)

        if NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=outfilename, auto_open=_AUTOOPEN)
            logging.info("Output written on {}".format(outfilename))

    def conditional_hist(self, residue, atom, atomlist, cslist, filtered=True, sd_limit=10, normalized=False,
                         ambiguity='*', outfilename='hist.html'):
        """
        Generates chemical shift histogram, which depends on the chemical shift values of given list of atoms
        in the same amino acid
        :param residue: 3 letter amino acid code
        :param atom: IUPAC atom name
        :param atomlist: known atom list as list
        :param cslist: corresponding chemical shift list as list
        :param filtered: True/False Filters based on standard deviation cutoff Default:True
        :param sd_limit: Number of time Standard deviation for filtering default: 10
        :param normalized: True/False Plots either Count/Density default: False
        :param outfilename: output file name
        :return: writes output in a html file
        """
        logging.info("This may take time to gather data from entire BMRB database!")
        if normalized:
            count = 'Density'
        else:
            count = 'Count'
        layout = plotly.graph_objs.Layout(
            barmode='overlay',
            xaxis=dict(autorange='reversed', title='Chemical Shift [ppm]'),
            yaxis=dict(title=count))
        data = [self.get_histogram_api(residue, atom, filtered, sd_limit, normalized,ambiguity)[0],
                self.get_conditional_histogram_api(residue, atom, atomlist, cslist, filtered, sd_limit, normalized,ambiguity)
                ]
        fig = plotly.graph_objs.Figure(data=data, layout=layout)

        if NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=outfilename, auto_open=_AUTOOPEN)
            logging.info("Output written on {}".format(outfilename))


def _called_directly():
    """ Figure out what to do if we were called on the command line
        rather than imported as a module."""

    # Specify some basic information about our command
    optparser = optparse.OptionParser(usage="usage: %prog",
                                      version=__version__,
                                      description="BMRB data visualization using python plotly")
    optparser.add_option("--hsqc", metavar="BMRBID", action="store",
                         dest="hsqc", default=None, type="string", nargs=1,
                         help="Outputs HSQC spectrum for the given BMRB entry.")
    optparser.add_option("--hist", metavar="residue atom", action="store",
                         dest="hist", default=None, type="string", nargs=2,
                         help="Outputs Histogram")
    optparser.add_option("--seq", metavar="One letter sequence", action="store",
                         dest="seq", default=None, nargs=1, type="string",
                         help="Random coil HSQC")
    optparser.add_option("--out", metavar="filename", action="store",
                         dest="outfile", default=None, nargs=1, type="string",
                         help="Output filename")

    # Options, parse 'em
    (options, cmd_input) = optparser.parse_args()

    if len(cmd_input) > 0:
        print("No arguments are allowed. Please see the options using --help.")
        sys.exit(0)
    # Check for command misuse
    if sum(1 for x in [options.hsqc,
                       options.hist, options.seq] if x) != 1:
        print("You have the following options \n"
              "python csviz.py --hsqc <BMRBID> --out <output file name>\n"
              "python csviz.py --hist <residue> <atom> --out <output file name>\n"
              "python csviz.py --seq <one letter sequence> --out  <output filename>\n"
              "python csviz.py --help")
        sys.exit(1)
    if options.outfile is None:
        print("Output file name must be specified with --out")
        sys.exit(1)
    if options.hsqc is not None:
        print(options.hsqc)
        s = Spectra()
        s.n15hsqc(bmrbid=options.hsqc.split(','), outfilename=options.outfile)
    elif options.hist is not None:
        h = Histogram()
        h.single_atom(residue=options.hist[0], atom=options.hist[1], outfilename=options.outfile)
    elif options.seq is not None:
        s = Spectra()
        s.n15hsqc(seq=options.seq, outfilename=options.outfile)
    else:
        print("Nothing specified")
        sys.exit(0)


if __name__ == "__main__":
    _called_directly()

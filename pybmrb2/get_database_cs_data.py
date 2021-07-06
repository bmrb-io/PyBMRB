from __future__ import print_function
import logging
import os
import sys
import ntpath
import json
import numpy

PY3 = (sys.version_info[0] == 3)
(scriptPath, scriptName) = ntpath.split(os.path.realpath(__file__))

if PY3:
    from urllib.request import urlopen, Request
else:
    from urllib2 import urlopen, Request

# Set the log level to INFO
#logging.getLogger().setLevel(logging.INFO)

_API_URL = "http://dev-api.bmrb.io/v2"


three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                     'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                     'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                     'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])

class ChemicalShiftStatistics(object):
    def __init__(self):
        pass

    @classmethod
    def get_dta_from_api(self,residue, atom):
        if residue == "*" and atom == "*":
            logging.error("Getting full database will overload the memory! Please chose a residue or atom.")
            raise ValueError("Getting full database will overload the memory! Please chose a residue or atom.")
        elif residue == "*":
            url = Request(_API_URL + "/search/chemical_shifts?atom_id={}".format(atom))
        elif atom == "*":
            url = Request(_API_URL + "/search/chemical_shifts?comp_id={}".format(residue))
        else:
            url = Request(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom))
        print (_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom))
        url.add_header('Application', 'PyBMRB')
        r = urlopen(url)
        dump = json.loads(r.read())
        return dump

    @classmethod
    def get_data(self,residue,atom,filtered=True, sd_limit=10, normalized=False, ambiguity='*',ph_min=None,ph_max=None,t_min=None,t_max=None):
        api_data=self.get_dta_from_api(residue=residue,atom=atom)
        cs_index = api_data['columns'].index('Atom_chem_shift.Val')
        ph_index = api_data['columns'].index('Sample_conditions.pH')
        temp_index = api_data['columns'].index('Sample_conditions.Temperature_K')
        res_index = api_data['columns'].index('Atom_chem_shift.Comp_ID')
        atm_index = api_data['columns'].index('Atom_chem_shift.Atom_ID')
        amb_index = api_data['columns'].index('Atom_chem_shift.Ambiguity_code')
        data=api_data['data']
        if ambiguity != '*':
            data=[i for i in data if i[amb_index]==ambiguity]
        if filtered:
            cs = [i[cs_index] for i in data]
            sd=numpy.std(cs)
            data=[i for i in data if -(sd_limit*sd)<i[cs_index]<(sd_limit*sd)]
        if ph_min is not None:
            data = [i for i in data if i[ph_index] >= ph_min]
        if ph_max is not None:
            data = [i for i in data if i[ph_index] <= ph_max]
        if t_min is not None:
            data = [i for i in data if i[temp_index] >= t_min]
        if t_max is not None:
            data = [i for i in data if i[temp_index] <= t_max]
        def transpose_list_of_lists(x):
            out_list=[[] for i in range(len(x[0]))]
            for i in x:
                for j in range(len(i)):
                    out_list[j].append(i[j])
            return out_list
        dat=transpose_list_of_lists(data)
        cs_data = [dat[cs_index],dat[ph_index],dat[temp_index],dat[res_index],dat[atm_index]]
        return cs_data





if __name__=="__main__":
    ChemicalShiftStatistics.get_data('ALA','HA',filtered=True,ambiguity=1)
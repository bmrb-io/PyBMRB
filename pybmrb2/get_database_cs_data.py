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
logging.getLogger().setLevel(logging.INFO)

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
        logging.info('Fetching chemical shift data for {}-{}'.format(residue,atom))
        url.add_header('Application', 'PyBMRB')
        r = urlopen(url)
        dump = json.loads(r.read())
        return dump

    @classmethod
    def get_data(self,residue,atom,filtered=True, sd_limit=10, ambiguity='*',
                 ph_min=None,ph_max=None,t_min=None,t_max=None,standard_amino_acids=True):
        standard = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS',
                    'ASP', 'SER', 'LYS', 'PRO', 'ASN',
                    'VAL', 'THR', 'HIS', 'TRP', 'PHE',
                    'ALA', 'MET', 'LEU', 'ARG', 'TYR',
                    'A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
        api_data=self.get_dta_from_api(residue=residue,atom=atom)
        columns = api_data['columns']
        cs_index = columns.index('Atom_chem_shift.Val')
        ph_index = columns.index('Sample_conditions.pH')
        temp_index = columns.index('Sample_conditions.Temperature_K')
        res_index = columns.index('Atom_chem_shift.Comp_ID')
        atm_index = columns.index('Atom_chem_shift.Atom_ID')
        amb_index = columns.index('Atom_chem_shift.Ambiguity_code')

        data=api_data['data']
        if len(data):
            if standard_amino_acids:
                data = [i for i in data if i[res_index] in standard]
            if ambiguity != '*':
                data=[i for i in data if i[amb_index]==ambiguity]
            if filtered:
                cs = [i[cs_index] for i in data]
                sd=numpy.std(cs)
                avg = numpy.mean(cs)
                min_cs = avg - (sd_limit*sd)
                max_cs = avg + (sd_limit*sd)
                data=[i for i in data if min_cs<i[cs_index]<max_cs]
            if ph_min is not None:
                data = [i for i in data if i[ph_index] is not None and i[ph_index] >= ph_min]
            if ph_max is not None:
                data = [i for i in data if i[ph_index] is not None and i[ph_index] <= ph_max]
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


        else:
            logging.info('Data for {}-{} not found in BMRB database'.format(residue,atom))
        #dat=transpose_list_of_lists(data)
        #cs_data = [dat[cs_index],dat[ph_index],dat[temp_index],dat[res_index],dat[atm_index]]
        return columns, data

    @classmethod
    def get_data_from_bmrb(self,residue=None,atom=None,list_of_atoms=None,filtered=True, sd_limit=10, ambiguity='*',
                 ph_min=None,ph_max=None,t_min=None,t_max=None,standard_amino_acids=True):
        if residue is None and atom is None and list_of_atoms is None:
            logging.error('Please provide residue name or atom name or list of atoms as a list ')
            raise TypeError('Please provide at leaset one of the three positional arguments: residue (or) atom (or) list_of_atoms')
        columns=[]
        out_dat=[]
        if list_of_atoms is not None:
            if type(list_of_atoms) is list:
                for atoms in list_of_atoms:
                    res = atoms.split("-")[0]
                    atm = atoms.split("-")[1]
                    cs_dat = self.get_data(residue=res,atom=atm,filtered=filtered,
                                           sd_limit=sd_limit,
                                           ambiguity=ambiguity,ph_min=ph_min,ph_max=ph_max,
                                           t_min=t_min,t_max=t_max,standard_amino_acids=standard_amino_acids)
                    if len(columns):
                        if cs_dat[0] == columns:
                            out_dat+=cs_dat[1]
                        else:
                            raise ValueError('BMRB API gives different column values')
                    else:
                        columns = cs_dat[0]
                        out_dat+=cs_dat[1]
            else:
                res = list_of_atoms.split("-")[0]
                atm = list_of_atoms.split("-")[1]
                cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                       sd_limit=sd_limit,
                                       ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                       t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat+=cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat+=cs_dat[1]
        else:
            logging.info("list of atoms was not provided")

        if type(residue) is list and type(atom) is list:
            for res in residue:
                for atm in atom:
                    cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                           sd_limit=sd_limit,
                                           ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                           t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
                    if len(columns):
                        if cs_dat[0] == columns:
                            out_dat+=cs_dat[1]
                        else:
                            raise ValueError('BMRB API gives different column values')
                    else:
                        columns = cs_dat[0]
                        out_dat+=cs_dat[1]
        elif type(residue) is list and atom is None:
            for res in residue:
                atm = '*'
                cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                       sd_limit=sd_limit,
                                       ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                       t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat+=cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat+=cs_dat[1]
        elif type(residue) is list and atom is not None:
            for res in residue:
                atm = atom
                cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                       sd_limit=sd_limit,
                                       ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                       t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat+=cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat+=cs_dat[1]
        elif residue is None and type(atom) is list:
            for atm in atom:
                res = '*'
                cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                       sd_limit=sd_limit,
                                       ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                       t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat+=cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat+=cs_dat[1]
        elif residue is not None and type(atom) is list:
            for atm in atom:
                res = residue
                cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                       sd_limit=sd_limit,
                                       ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                       t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat+=cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat+=cs_dat[1]
        elif residue is not None and atom is not None:
            res = residue
            atm = atom
            cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                   sd_limit=sd_limit,
                                   ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                   t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
            if len(columns):
                if cs_dat[0] == columns:
                    out_dat+=cs_dat[1]
                else:
                    raise ValueError('BMRB API gives different column values')
            else:
                columns = cs_dat[0]
                out_dat+=cs_dat[1]
        elif residue is None and atom is not None:
            res = '*'
            atm = atom
            cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                   sd_limit=sd_limit,
                                   ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                   t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
            if len(columns):
                if cs_dat[0] == columns:
                    out_dat+=cs_dat[1]
                else:
                    raise ValueError('BMRB API gives different column values')
            else:
                columns = cs_dat[0]
                out_dat+=cs_dat[1]
        elif residue is not None and atom is None:
            res = residue
            atm = '*'
            cs_dat = self.get_data(residue=res, atom=atm, filtered=filtered,
                                   sd_limit=sd_limit,
                                   ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                   t_min=t_min, t_max=t_max,standard_amino_acids=standard_amino_acids)
            if len(columns):
                if cs_dat[0] == columns:
                    out_dat+=cs_dat[1]
                else:
                    raise ValueError('BMRB API gives different column values')
            else:
                columns = cs_dat[0]
                out_dat+=cs_dat[1]
        else:
            logging.info('residue and atom values were not provided')
        return columns, out_dat
    @classmethod
    def list_do_dict(self,columns,data):
        entity_index=columns.index('Atom_chem_shift.Entity_ID')
        entry_index = columns.index('Atom_chem_shift.Entry_ID')
        cs_index = columns.index('Atom_chem_shift.Val')
        res_index = columns.index('Atom_chem_shift.Comp_ID')
        seq_index = columns.index('Atom_chem_shift.Comp_index_ID')
        atm_index = columns.index('Atom_chem_shift.Atom_ID')
        list_index = columns.index('Atom_chem_shift.Assigned_chem_shift_list_ID')
        amb_index = columns.index('Atom_chem_shift.Ambiguity_code')
        cs_dict={}
        for dat in data:
            entry_id = dat[entry_index]
            entity_id = dat[entity_index]
            res = dat[res_index]
            seq = dat[seq_index]
            list_id = dat[list_index]
            key = '{}-{}-{}-{}-{}'.format(entry_id,entity_id,seq,res,list_id)
            if key not in cs_dict.keys():
                cs_dict[key]={}
            atom = dat[atm_index]
            cs = dat[cs_index]
            ambi = dat[amb_index]
            if atom not in cs_dict[key].keys():
                cs_dict[key][atom]=(cs,ambi)
            else:
                logging.warning('Duplicate key found {}'.format(dat))
        return cs_dict

    @classmethod
    def get_2d_chemical_shifts(self,residue,atom1,atom2,filtered=True, sd_limit=10,
                               ambiguity1='*', ambiguity2='*',
                               ph_min=None,ph_max=None,t_min=None,t_max=None,standard_amino_acids=True):
        x=[]
        y=[]
        cs_data = self.get_data_from_bmrb(residue=residue,ph_min=ph_min,ph_max=ph_max,
                                          t_min=t_min,t_max=t_max,standard_amino_acids=standard_amino_acids)
        cs_dict = self.list_do_dict(cs_data[0],cs_data[1])
        for key in cs_dict.keys():
            if atom1 in cs_dict[key].keys() and atom2 in cs_dict[key].keys():
                if ambiguity1 !='*' or ambiguity2 !='*':
                    if ambiguity1 != '*':
                        if cs_dict[key][atom2][1] == ambiguity2:
                            x.append(cs_dict[key][atom1][0])
                            y.append(cs_dict[key][atom2][0])
                    elif ambiguity2 != '*':
                        if cs_dict[key][atom1][1] == ambiguity1:
                            x.append(cs_dict[key][atom1][0])
                            y.append(cs_dict[key][atom2][0])
                    else:
                        if cs_dict[key][atom1][1] == ambiguity1 and cs_dict[key][atom2][1] == ambiguity2:
                            x.append(cs_dict[key][atom1][0])
                            y.append(cs_dict[key][atom2][0])
                else:
                    x.append(cs_dict[key][atom1][0])
                    y.append(cs_dict[key][atom2][0])
        cs_x = []
        cs_y = []
        if filtered:
            mean_x = numpy.mean(x)
            mean_y = numpy.mean(y)
            std_x = numpy.std(x)
            std_y = numpy.std(y)
            min_x = mean_x - std_x*sd_limit
            max_x = mean_x + std_x*sd_limit
            min_y = mean_y - std_y*sd_limit
            max_y = mean_y + std_y*sd_limit
            for i in range(len(x)):
                if min_x<=x[i]<=max_x and min_y<=y[i]<=max_y:
                    cs_x.append(x[i])
                    cs_y.append(y[i])
        return cs_x,cs_y

    @classmethod
    def get_filtered_data_from_bmrb(self, residue, atom, filtering_rules,
                           ph_min=None, ph_max=None, t_min=None, t_max=None, standard_amino_acids=True):
        cs_data = self.get_data_from_bmrb(residue=residue, ph_min=ph_min, ph_max=ph_max,
                                          t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
        def filter(cs_dict,atm,cs_val):
            if 'H' in atm:
                cs_width = 0.1
            if 'C' in atm:
                cs_width = 2.0
            if 'N' in atm:
                cs_width = 2.0
            out_cs_dict = {}
            for key in cs_dict.keys():
                try:
                    if cs_val-cs_width <= cs_dict[key][atm][0] <= cs_val+cs_width:
                        out_cs_dict[key] = cs_dict[key]
                except KeyError:
                    pass
            return out_cs_dict
        cs_dict = self.list_do_dict(cs_data[0], cs_data[1])
        for rule in filtering_rules:
            cs_dict = filter(cs_dict,rule[0],rule[1])
        x=[]
        for key in cs_dict.keys():
            for atm in cs_dict[key].keys():
                if atm == atom:
                    x.append(cs_dict[key][atm][0])
        return x















if __name__=="__main__":
    # x=ChemicalShiftStatistics.get_data_from_bmrb(list_of_atoms='ALA-N')
    # #print (x[0])
    # y=ChemicalShiftStatistics.list_do_dict(x[0],x[1])
    #x=ChemicalShiftStatistics.get_2d_chemical_shifts(residue='ALA',atom1='CA',atom2='CB')
    x=ChemicalShiftStatistics.get_filtered_data_from_bmrb(residue='THR',atom='N',filtering_rules=[('CB',69.51),('CA',60.79),('H',8.13)])
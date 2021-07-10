#!/usr/bin/env python3
import json
import logging
from urllib.request import urlopen, Request
import numpy
# Set the log level to INFO
logging.getLogger().setLevel(logging.INFO)

# _API_URL = "http://dev-api.bmrb.io/v2"
_API_URL = "http://api.bmrb.io/v2"

three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                     'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                     'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                     'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])


class ChemicalShiftStatistics(object):
    '''
    Fetches chemical shift data from BMRB using BMRB-API
    '''

    def __init__(self):
        pass

    @staticmethod
    def _get_data_from_api(residue, atom):
        '''
        Dumps the BMRB-API return for given residue and atom. Not intend to call directly

        :param residue: Residue name in IUPAC format
        :param atom: Atom name in IUPAC format
        :return: API dump
        '''
        # logging.info('Fetching chemical shift data for {}-{}'.format(residue,atom))
        if residue == "*" and atom == "*":
            logging.error("Getting full database will overload the memory! Please chose a residue or atom.")
            raise ValueError("Getting full database will overload the memory! Please chose a residue or atom.")
        elif residue == "*":
            url = Request(_API_URL + "/search/chemical_shifts?atom_id={}".format(atom))
        elif atom == "*":
            url = Request(_API_URL + "/search/chemical_shifts?comp_id={}".format(residue))
        else:
            url = Request(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom))
        logging.debug(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom))
        logging.info('Fetching chemical shift data for {}-{}'.format(residue, atom))
        url.add_header('Application', 'PyBMRB')
        r = urlopen(url)
        dump = json.loads(r.read())
        return dump

    @classmethod
    def get_data(cls, residue, atom, filtered=True, sd_limit=10, ambiguity='*',
                 ph_min=None, ph_max=None, t_min=None, t_max=None, standard_amino_acids=True):
        '''
        Fetches the data from BMRB-API for given residue and atom and filters based of the values provided in the parameters

        :param residue: Residue name in IPUPAC format; use '*' for any residue
        :param atom: Atom name in IUPAC format; Wild card supported '*' for any atom
        :param filtered: Filters values beyond (sd_limt)*(standard deviation) on both sides of the mean; default:True
        :param sd_limit: scaling factor used to filter data based on standard deviation; default 10
        :param ambiguity: ambiguity filter; default '*' => no filter
        :param ph_min: PH filter (min);default None
        :param ph_max: PH filter (max); default None
        :param t_min: Temperature filter (min); default None
        :param t_max: Temperature filter (max); default None
        :param standard_amino_acids: get data only form 20 natural amino acids,4 standard DNA and 4 standard RNA; default:True
        :return: column names and data as tuple (columns,data)
        '''
        standard = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS',
                    'ASP', 'SER', 'LYS', 'PRO', 'ASN',
                    'VAL', 'THR', 'HIS', 'TRP', 'PHE',
                    'ALA', 'MET', 'LEU', 'ARG', 'TYR',
                    'A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
        api_data = cls._get_data_from_api(residue=residue, atom=atom)
        columns = api_data['columns']
        cs_index = columns.index('Atom_chem_shift.Val')
        ph_index = columns.index('Sample_conditions.pH')
        temp_index = columns.index('Sample_conditions.Temperature_K')
        res_index = columns.index('Atom_chem_shift.Comp_ID')
        atm_index = columns.index('Atom_chem_shift.Atom_ID')
        amb_index = columns.index('Atom_chem_shift.Ambiguity_code')
        data = api_data['data']
        if len(data):
            if standard_amino_acids:
                data = [i for i in data if i[res_index] in standard]
            if ambiguity != '*':
                data = [i for i in data if i[amb_index] == ambiguity]
            if filtered:
                cs = [i[cs_index] for i in data]
                sd = numpy.std(cs)
                avg = numpy.mean(cs)
                min_cs = avg - (sd_limit * sd)
                max_cs = avg + (sd_limit * sd)
                data = [i for i in data if min_cs < i[cs_index] < max_cs]
            if ph_min is not None:
                data = [i for i in data if i[ph_index] is not None and i[ph_index] >= ph_min]
            if ph_max is not None:
                data = [i for i in data if i[ph_index] is not None and i[ph_index] <= ph_max]
            if t_min is not None:
                data = [i for i in data if i[temp_index] >= t_min]
            if t_max is not None:
                data = [i for i in data if i[temp_index] <= t_max]
        else:
            logging.info('Data for {}-{} not found in BMRB database'.format(residue, atom))
        return columns, data

    @classmethod
    def get_data_from_bmrb(cls, residue=None, atom=None, list_of_atoms=None, filtered=True, sd_limit=10, ambiguity='*',
                           ph_min=None, ph_max=None, t_min=None, t_max=None, standard_amino_acids=True):
        '''
        Fetches the data from BMRB-API for given residue/list of residues (and/or)  atom/list of atoms and
         filters based of the values provided in the parameters

        :param residue: Single residue name or list of residue names  in IUPAC format
        :param atom: Single atom name or list of atom names in IUPAC format
        :param list_of_atoms: list of atoms; example ['ALA-CB','CYS-N','TYR-CB']
        :param filtered: Filters values beyond (sd_limt)*(standard deviation) on both sides of the mean; default:True
        :param sd_limit: scaling factor used to filter data based on standard deviation; default 10
        :param ambiguity: ambiguity filter; default '*' => no filter
        :param ph_min: PH filter (min);default None
        :param ph_max: PH filter (max); default None
        :param t_min: Temperature filter (min); default None
        :param t_max: Temperature filter (max); default None
        :param standard_amino_acids: get data only form 20 natural amino acids,4 standard DNA and 4 standard RNA; default:True
        :return: column names and data as tuple (columns,data)
        '''
        if residue is None and atom is None and list_of_atoms is None:
            logging.error('Please provide residue name or atom name or list of atoms as a list ')
            raise TypeError(
                'Please provide at leaset one of the three positional arguments: residue (or) atom (or) list_of_atoms')
        columns = []
        out_dat = []
        if list_of_atoms is not None:
            if type(list_of_atoms) is list:
                for atoms in list_of_atoms:
                    res = atoms.split("-")[0]
                    atm = atoms.split("-")[1]
                    cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                          sd_limit=sd_limit,
                                          ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                          t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
                    if len(columns):
                        if cs_dat[0] == columns:
                            out_dat += cs_dat[1]
                        else:
                            raise ValueError('BMRB API gives different column values')
                    else:
                        columns = cs_dat[0]
                        out_dat += cs_dat[1]
            else:
                res = list_of_atoms.split("-")[0]
                atm = list_of_atoms.split("-")[1]
                cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                      sd_limit=sd_limit,
                                      ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                      t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat += cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat += cs_dat[1]
        else:
            logging.debug("list of atoms was not provided")

        if type(residue) is list and type(atom) is list:
            for res in residue:
                for atm in atom:
                    cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                          sd_limit=sd_limit,
                                          ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                          t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
                    if len(columns):
                        if cs_dat[0] == columns:
                            out_dat += cs_dat[1]
                        else:
                            raise ValueError('BMRB API gives different column values')
                    else:
                        columns = cs_dat[0]
                        out_dat += cs_dat[1]
        elif type(residue) is list and atom is None:
            for res in residue:
                atm = '*'
                cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                      sd_limit=sd_limit,
                                      ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                      t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat += cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat += cs_dat[1]
        elif type(residue) is list and atom is not None:
            for res in residue:
                atm = atom
                cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                      sd_limit=sd_limit,
                                      ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                      t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat += cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat += cs_dat[1]
        elif residue is None and type(atom) is list:
            for atm in atom:
                res = '*'
                cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                      sd_limit=sd_limit,
                                      ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                      t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat += cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat += cs_dat[1]
        elif residue is not None and type(atom) is list:
            for atm in atom:
                res = residue
                cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                      sd_limit=sd_limit,
                                      ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                      t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
                if len(columns):
                    if cs_dat[0] == columns:
                        out_dat += cs_dat[1]
                    else:
                        raise ValueError('BMRB API gives different column values')
                else:
                    columns = cs_dat[0]
                    out_dat += cs_dat[1]
        elif residue is not None and atom is not None:
            res = residue
            atm = atom
            cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                  sd_limit=sd_limit,
                                  ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                  t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
            if len(columns):
                if cs_dat[0] == columns:
                    out_dat += cs_dat[1]
                else:
                    raise ValueError('BMRB API gives different column values')
            else:
                columns = cs_dat[0]
                out_dat += cs_dat[1]
        elif residue is None and atom is not None:
            res = '*'
            atm = atom
            cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                  sd_limit=sd_limit,
                                  ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                  t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
            if len(columns):
                if cs_dat[0] == columns:
                    out_dat += cs_dat[1]
                else:
                    raise ValueError('BMRB API gives different column values')
            else:
                columns = cs_dat[0]
                out_dat += cs_dat[1]
        elif residue is not None and atom is None:
            res = residue
            atm = '*'
            cs_dat = cls.get_data(residue=res, atom=atm, filtered=filtered,
                                  sd_limit=sd_limit,
                                  ambiguity=ambiguity, ph_min=ph_min, ph_max=ph_max,
                                  t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)
            if len(columns):
                if cs_dat[0] == columns:
                    out_dat += cs_dat[1]
                else:
                    raise ValueError('BMRB API gives different column values')
            else:
                columns = cs_dat[0]
                out_dat += cs_dat[1]
        else:
            logging.info('residue and atom values were not provided')
        return columns, out_dat

    @staticmethod
    def _list_to_dict(columns, data):
        '''
        Converts the outputs from get_data_from_bmrb into a dictionary. Not intend to call directly

        :param columns: Column headers
        :param data: data as list of lists
        :return: chemical shift dictionary; { entry_id-entity_id-seq_id-residue-chemical_shift_list_id:{atom:(chemical shift, ambiguity code}}
        '''
        entity_index = columns.index('Atom_chem_shift.Entity_ID')
        entity_assembly_index = columns.index('Atom_chem_shift.Entity_assembly_ID')
        entry_index = columns.index('Atom_chem_shift.Entry_ID')
        cs_index = columns.index('Atom_chem_shift.Val')
        res_index = columns.index('Atom_chem_shift.Comp_ID')
        seq_index = columns.index('Atom_chem_shift.Comp_index_ID')
        atm_index = columns.index('Atom_chem_shift.Atom_ID')
        list_index = columns.index('Atom_chem_shift.Assigned_chem_shift_list_ID')
        amb_index = columns.index('Atom_chem_shift.Ambiguity_code')
        cs_dict = {}
        for dat in data:
            entry_id = dat[entry_index]
            entity_id = dat[entity_index]
            entity_assembly_id = dat[entity_assembly_index]
            res = dat[res_index]
            seq = dat[seq_index]
            list_id = dat[list_index]
            key = '{}-{}-{}-{}-{}'.format(entry_id, entity_assembly_id, seq, res, list_id)
            if key not in cs_dict.keys():
                cs_dict[key] = {}
            atom = dat[atm_index]
            cs = dat[cs_index]
            ambi = dat[amb_index]
            if atom not in cs_dict[key].keys():
                cs_dict[key][atom] = (cs, ambi)
            else:
                logging.warning('Duplicate key found {}'.format(dat))
        return cs_dict

    @classmethod
    def get_2d_chemical_shifts(cls, residue, atom1, atom2, filtered=True, sd_limit=10,
                               ambiguity1='*', ambiguity2='*',
                               ph_min=None, ph_max=None, t_min=None, t_max=None):
        '''
        Fetches chemical shift data for a given residue and combines the desired two atoms chemical shift value from the
        same residue as two dimensional list

        :param residue: residue name in IUPAC format
        :param atom1: atom name in IUPAC format
        :param atom2: atom name in IPUPAC format
        :param filtered: Filters values beyond (sd_limt)*(standard deviation) on both sides of the mean; default:True
        :param sd_limit: scaling factor used to filter data based on standard deviation; default 10
        :param ambiguity1: ambiguity filter; default '*' (no filter)
        :param ambiguity2: ambiguity filter; default '*'  (no filter)
        :param ph_min: PH filter (min);default None
        :param ph_max: PH filter (max); default None
        :param t_min: Temperature filter (min); default None
        :param t_max: Temperature filter (max); default None
        :return: tuple of lists (atom1_cs,atom2_cs as list)
        '''
        x = []
        y = []
        cs_data = cls.get_data_from_bmrb(residue=residue, ph_min=ph_min, ph_max=ph_max,
                                         t_min=t_min, t_max=t_max, standard_amino_acids=False)
        cs_dict = cls._list_to_dict(cs_data[0], cs_data[1])
        for key in cs_dict.keys():
            if atom1 in cs_dict[key].keys() and atom2 in cs_dict[key].keys():
                if ambiguity1 != '*' or ambiguity2 != '*':
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
            min_x = mean_x - std_x * sd_limit
            max_x = mean_x + std_x * sd_limit
            min_y = mean_y - std_y * sd_limit
            max_y = mean_y + std_y * sd_limit
            for i in range(len(x)):
                if min_x <= x[i] <= max_x and min_y <= y[i] <= max_y:
                    cs_x.append(x[i])
                    cs_y.append(y[i])
        return cs_x, cs_y

    @classmethod
    def get_filtered_data_from_bmrb(cls, residue, atom, filtering_rules,
                                    ph_min=None, ph_max=None, t_min=None, t_max=None, standard_amino_acids=True):
        '''
        Fetches the chemical shift data for a given residue and filters them using filtering rules

        :param residue: residue name in IUPAC format; example 'CYS'
        :param atom: atom name in IUPAC format; example 'CB'
        :param filtering_rules: list of atoms and chemical shift values as tuples; example[('CA',64.5),('H',7.8)]
        :param ph_min: PH filter (min);default None
        :param ph_max: PH filter (max); default None
        :param t_min: Temperature filter (min); default None
        :param t_max: Temperature filter (max); default None
        :param standard_amino_acids: get data only form 20 natural amino acids,4 standard DNA and 4 standard RNA; default:True
        :return: chemical shift values as list
        '''
        cs_data = cls.get_data_from_bmrb(residue=residue, ph_min=ph_min, ph_max=ph_max,
                                         t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)

        def filter(cs_dict, atm, cs_val):
            if 'H' in atm:
                cs_width = 0.1
            if 'C' in atm:
                cs_width = 2.0
            if 'N' in atm:
                cs_width = 2.0
            out_cs_dict = {}
            for key in cs_dict.keys():
                try:
                    if cs_val - cs_width <= cs_dict[key][atm][0] <= cs_val + cs_width:
                        out_cs_dict[key] = cs_dict[key]
                except KeyError:
                    pass
            return out_cs_dict

        cs_dict = cls._list_to_dict(cs_data[0], cs_data[1])
        for rule in filtering_rules:
            cs_dict = filter(cs_dict, rule[0], rule[1])
        x = []
        for key in cs_dict.keys():
            for atm in cs_dict[key].keys():
                if atm == atom:
                    x.append(cs_dict[key][atm][0])
        return x

    @classmethod
    def get_statistics(cls, residue=None, atom=None, list_of_atoms=None, filtered=True, sd_limit=10, ambiguity='*',
                       ph_min=None, ph_max=None, t_min=None, t_max=None, standard_amino_acids=True, verbose=False):
        '''
        Provides chemical shift statistics like Mean, Median, Standard deviation for a given residue / atom /list of atoms

        :param residue: Single residue name or list of residue names  in IUPAC format
        :param atom: Single atom name or list of atom names in IUPAC format
        :param list_of_atoms: list of atoms; example ['ALA-CB','CYS-N','TYR-CB']
        :param filtered: Filters values beyond (sd_limt)*(standard deviation) on both sides of the mean; default:True
        :param sd_limit: scaling factor used to filter data based on standard deviation; default 10
        :param ambiguity: ambiguity filter; default '*' => no filter
        :param ph_min: PH filter (min);default None
        :param ph_max: PH filter (max); default None
        :param t_min: Temperature filter (min); default None
        :param t_max: Temperature filter (max); default None
        :param standard_amino_acids: get data only form 20 natural amino acids,4 standard DNA and 4 standard RNA; default:True
        :param verbose: display the statistics on screen; default False
        :return: chemcial shift statistics as a dictionary
        '''
        if residue is None and atom is None and list_of_atoms is None:
            logging.error('Please provide residue name or atom name or list of atoms as a list ')
            raise TypeError(
                'Please provide at leaset one of the three positional arguments: residue (or) atom (or) list_of_atoms')
        columns, data = cls.get_data_from_bmrb(residue=residue,
                                               atom=atom,
                                               list_of_atoms=list_of_atoms,
                                               filtered=filtered,
                                               sd_limit=sd_limit,
                                               ambiguity=ambiguity,
                                               ph_min=ph_min,
                                               ph_max=ph_max,
                                               t_min=t_min,
                                               t_max=t_max,
                                               standard_amino_acids=standard_amino_acids)
        res_index = columns.index('Atom_chem_shift.Comp_ID')
        atm_index = columns.index('Atom_chem_shift.Atom_ID')
        cs_index = columns.index('Atom_chem_shift.Val')
        atm_list = list(set(['{}-{}'.format(i[res_index], i[atm_index]) for i in data]))
        cs_stat = {}
        for atm in atm_list:
            cs_stat[atm] = {}
            res = atm.split("-")[0]
            at = atm.split("-")[1]
            x = [i[cs_index] for i in data if i[res_index] == res and i[atm_index] == at]
            cs_stat[atm]['cs'] = x
            cs_stat[atm]['mean'] = round(numpy.mean(x), 3)
            cs_stat[atm]['median'] = round(numpy.median(x), 3)
            cs_stat[atm]['std'] = round(numpy.std(x), 3)
            cs_stat[atm]['min'] = min(x)
            cs_stat[atm]['max'] = max(x)
            cs_stat[atm]['count'] = len(x)
        if verbose:
            for key in cs_stat.keys():
                print('\t\t{}'.format(key))
                print('\t\tMean\t\t\t\t:{}'.format(cs_stat[key]['mean']))
                print('\t\tMedian\t\t\t\t:{}'.format(cs_stat[key]['median']))
                print('\t\tStandard deviation\t:{}'.format(cs_stat[key]['std']))
                print('\t\tMinimum\t\t\t\t:{}'.format(cs_stat[key]['min']))
                print('\t\tMaximum\t\t\t\t:{}'.format(cs_stat[key]['max']))
                print('\t\tCount\t\t\t\t:{}'.format(cs_stat[key]['count']))
                print("\n")
        return cs_stat

# if __name__=="__main__":
#     # x=ChemicalShiftStatistics.get_data_from_bmrb(list_of_atoms='ALA-N')
#     # #print (x[0])
#     # y=ChemicalShiftStatistics._list_to_dict(x[0],x[1])
#     #x=ChemicalShiftStatistics.get_2d_chemical_shifts(residue='ALA',atom1='CA',atom2='CB')
#     #x=ChemicalShiftStatistics.get_filtered_data_from_bmrb(residue='THR',atom='N',filtering_rules=[('CB',69.51),('CA',60.79),('H',8.13)])
#     ChemicalShiftStatistics.get_statistics(residue='GLY',atom='XX')

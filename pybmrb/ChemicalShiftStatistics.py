#!/usr/bin/env python3
"""
This module fetches chunks of chemical shift data from BMRB database for given list of atoms or residues.
It exports the chemical shift data for Histogram module. This module is not intended to call directly, but can be
used for data mining.
"""
import json
import logging
from urllib.request import urlopen, Request
import numpy
from typing import Union, List, Optional

# Set the log level to INFO
logging.getLogger().setLevel(logging.INFO)

# _API_URL = "http://dev-api.bmrb.io/v2"
_API_URL = "http://api.bmrb.io/v2"

three_letter_code = {'I': 'ILE', 'Q': 'GLN', 'G': 'GLY', 'E': 'GLU', 'C': 'CYS',
                     'D': 'ASP', 'S': 'SER', 'K': 'LYS', 'P': 'PRO', 'N': 'ASN',
                     'V': 'VAL', 'T': 'THR', 'H': 'HIS', 'W': 'TRP', 'F': 'PHE',
                     'A': 'ALA', 'M': 'MET', 'L': 'LEU', 'R': 'ARG', 'Y': 'TYR'}
one_letter_code = dict([(value, key) for key, value in three_letter_code.items()])


def _get_data_from_api(residue: str,
                       atom: str) -> dict:
    """
    Dumps the BMRB-API return for given residue and atom. Not intend to call directly

    :param residue: Residue name in IUPAC format
    :type residue: str
    :param atom: Atom name in IUPAC format
    :type atom: str
    :return: API dump
    :rtype: dict
    """
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
    logging.debug(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom))
    logging.info('Fetching chemical shift data for {}-{}'.format(residue, atom))
    url.add_header('Application', 'PyBMRB')
    r = urlopen(url)
    dump = json.loads(r.read())
    return dump


def get_data(residue: str,
             atom: str,
             filtered: Optional[bool] = True,
             sd_limit: Optional[float] = 10.0,
             ambiguity: Optional[Union[str, int]] = '*',
             ph_min: Optional[float] = None,
             ph_max: Optional[float] = None,
             t_min: Optional[float] = None,
             t_max: Optional[float] = None,
             standard_amino_acids: Optional[bool] = True) -> tuple:
    """
    Fetches the data from BMRB-API for given residue and atom and filters based of the values provided in the
    parameters

    :param residue: Residue name in IUPAC format; use \'*\' for any residue
    :type residue: str
    :param atom: Atom name in IUPAC format; Wild card supported; example HB*; use '*' for any atom
    :type atom: str
    :param filtered: Filters values beyond (sd_limit)*(standard deviation) on both sides of the mean,
        defaults to True
    :type filtered: bool, optional
    :param sd_limit: scaling factor used to filter data based on standard deviation, defaults to 10.0
    :type sd_limit: float, optional
    :param ambiguity: filter based on chemical shift ambiguity code, defaults to * means no filter
    :type ambiguity: str, optional
    :param ph_min: minimum value for the filter based on PH value, defaults to None
    :type ph_min: float, optional
    :param ph_max: maximum value for the filter based on PH value, defaults to None
    :type ph_max: float, optional
    :param t_min: minimum value for the filter based on temperature value, defaults to None
    :type t_min: float, optional
    :param t_max: maximum value for the filter based on temperature value, defaults to None
    :type t_max: float, optional
    :param standard_amino_acids: get data only form standard amino acids and nucleic acids when * is used for residue,
        defaults to True
    :type standard_amino_acids: bool, optional
    :return: column names and chemical shift data as tuple (columns,data)
    :rtype: tuple
    """

    standard = ['ILE', 'GLN', 'GLY', 'GLU', 'CYS',
                'ASP', 'SER', 'LYS', 'PRO', 'ASN',
                'VAL', 'THR', 'HIS', 'TRP', 'PHE',
                'ALA', 'MET', 'LEU', 'ARG', 'TYR',
                'A', 'C', 'G', 'U', 'DA', 'DC', 'DG', 'DT']
    api_data = _get_data_from_api(residue=residue, atom=atom)
    columns = api_data['columns']
    cs_index = columns.index('Atom_chem_shift.Val')
    ph_index = columns.index('Sample_conditions.pH')
    temp_index = columns.index('Sample_conditions.Temperature_K')
    res_index = columns.index('Atom_chem_shift.Comp_ID')
    atm_index = columns.index('Atom_chem_shift.Atom_ID')
    amb_index = columns.index('Atom_chem_shift.Ambiguity_code')
    data = api_data['data']
    for i in range(len(data)):
        data[i][cs_index] = float(data[i][cs_index])
        if t_min is not None or t_max is not None:
            try:
                data[i][temp_index] = float(data[i][temp_index])
            except TypeError:
                pass
        if ph_min is not None or ph_max is not None:
            try:
                data[i][ph_index] = float(data[i][ph_index])
            except TypeError:
                pass
    if len(data):
        if standard_amino_acids:
            data = [i for i in data if i[res_index] in standard]
        if ambiguity != '*':
            data = [i for i in data if i[amb_index] == ambiguity]
        if filtered:
            atom_list = set([i[atm_index] for i in data])
            a_dict = {}
            for atm in atom_list:
                a_dict[atm] = {}
                tmp_cs = [float(i[cs_index]) for i in data if i[atm_index] == atm]
                a_dict[atm]['sd'] = numpy.std(tmp_cs)
                a_dict[atm]['avg'] = numpy.mean(tmp_cs)
                a_dict[atm]['min_cs'] = a_dict[atm]['avg'] - (sd_limit * a_dict[atm]['sd'])
                a_dict[atm]['max_cs'] = a_dict[atm]['avg'] + (sd_limit * a_dict[atm]['sd'])

            # cs = [i[cs_index] for i in data]
            # sd = numpy.std(cs)
            # avg = numpy.mean(cs)
            # min_cs = avg - (sd_limit * sd)
            # max_cs = avg + (sd_limit * sd)
            data = [i for i in data if
                    a_dict[i[atm_index]]['min_cs'] < float(i[cs_index]) < a_dict[i[atm_index]]['max_cs']]
        if ph_min is not None:
            data2 = [i for i in data if i[ph_index] is not None]
            data = [i for i in data2 if i[ph_index] is not None and i[ph_index] >= ph_min]
        if ph_max is not None:
            data2 = [i for i in data if i[ph_index] is not None]
            data = [i for i in data2 if i[ph_index] is not None and i[ph_index] <= ph_max]
        if t_min is not None:
            data2 = [i for i in data if i[temp_index] is not None]
            data = [i for i in data2 if i[temp_index] >= t_min]
        if t_max is not None:
            data2 = [i for i in data if i[temp_index] is not None]
            data = [i for i in data2 if i[temp_index] <= t_max]
    else:
        logging.info('Data for {}-{} not found in BMRB database'.format(residue, atom))
    return columns, data


def get_data_from_bmrb(residue: str = None,
                       atom: str = None,
                       list_of_atoms: Optional[Union[str, List[str]]] = None,
                       filtered: Optional[bool] = True,
                       sd_limit: Optional[float] = 10.0,
                       ambiguity: Optional[Union[str, int]] = '*',
                       ph_min: Optional[float] = None,
                       ph_max: Optional[float] = None,
                       t_min: Optional[float] = None,
                       t_max: Optional[float] = None,
                       standard_amino_acids: Optional[bool] = True) -> tuple:
    """
    Fetches the data from BMRB-API for given residue/list of residues (and/or)  atom/list of atoms and
     filters based of the values provided in the parameters

    :param residue: Residue name in IUPAC format; use '*' for any residue
    :type residue: str
    :param atom: Atom name in IUPAC format; Wild card supported; example HB*; use '*' for any atom
    :type atom: str
    :param list_of_atoms: list of atoms; example ['ALA-CB','CYS-N','TYR-CB'], defaults to None
    :type list_of_atoms: list, optional
    :param filtered: Filters values beyond (sd_limit)*(standard deviation) on both sides of the mean,
        defaults to True
    :type filtered: bool, optional
    :param sd_limit: scaling factor used to filter data based on standard deviation, defaults to 10.0
    :type sd_limit: float, optional
    :param ambiguity: filter based on chemical shift ambiguity code, defaults to '*' means no filter
    :type ambiguity: str, optional
    :param ph_min: minimum value for the filter based on PH value, defaults to None
    :type ph_min: float, optional
    :param ph_max: maximum value for the filter based on PH value, defaults to None
    :type ph_max: float, optional
    :param t_min: minimum value for the filter based on temperature value, defaults to None
    :type t_min: float, optional
    :param t_max: maximum value for the filter based on temperature value, defaults to None
    :type t_max: float, optional
    :param standard_amino_acids: get data only form standard amino acids and nucleic acids when '*' is used for residue,
        defaults to True
    :type standard_amino_acids: bool, optional
    :return: column names and chemical shift data as tuple (columns,data)
    :rtype: tuple
    """

    if residue is None and atom is None and list_of_atoms is None:
        logging.error('Please provide residue name or atom name or list of atoms as a list ')
        raise TypeError(
            'Please provide at least one of the three positional arguments: residue (or) atom (or) list_of_atoms')
    columns = []
    out_dat = []
    if list_of_atoms is not None:
        if type(list_of_atoms) is list:
            for atoms in list_of_atoms:
                res = atoms.split("-")[0]
                atm = atoms.split("-")[1]
                cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
            cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
                cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
            cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
            cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
            cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
            cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
        cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
        cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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
        cs_dat = get_data(residue=res, atom=atm, filtered=filtered,
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


def _list_to_dict(columns: List[str],
                  data: list) -> dict:
    """
    Converts the outputs from get_data_from_bmrb into a dictionary. Not intended to call directly

    :param columns: Column headers
    :type columns: list
    :param data: data as list of lists
    :type data: list
    :return: chemical shift dictionary; { entry_id-entity_id-seq_id-residue-chemical_shift_list_id:
    {atom:(chemical shift, ambiguity code}}
    :rtype: dict
    """
    # entity_index = columns.index('Atom_chem_shift.Entity_ID')
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
        # entity_id = dat[entity_index]
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


def get_2d_chemical_shifts(residue: str,
                           atom1: str,
                           atom2: str,
                           filtered: Optional[bool] = True,
                           sd_limit: Optional[float] = 10,
                           ambiguity1: Optional[Union[str, int]] = '*',
                           ambiguity2: Optional[Union[str, int]] = '*',
                           ph_min: Optional[float] = None,
                           ph_max: Optional[float] = None,
                           t_min: Optional[float] = None,
                           t_max: Optional[float] = None) -> tuple:
    """
    Fetches chemical shift data for a given residue and combines the desired two atoms chemical shift value from the
    same residue as two dimensional list

    :param residue: residue name in IUPAC format
    :type residue: str
    :param atom1: atom name in IUPAC format
    :type atom1: str
    :param atom2: atom name in IUPAC format
    :type atom2: str
    :param filtered: Filters values beyond (sd_limit)*(standard deviation) on both sides of the mean,
        defaults to True
    :type filtered: bool, optional
    :param sd_limit: scaling factor used to filter data based on standard deviation, defaults to 10.0
    :type sd_limit: float, optional
    :param ambiguity1: filter based on chemical shift ambiguity code for atom1, defaults to '*' means no filter
    :type ambiguity1: str, optional
    :param ambiguity2: filter based on chemical shift ambiguity code for atom2, defaults to '*' means no filter
    :type ambiguity2: str, optional
    :param ph_min: minimum value for the filter based on PH value, defaults to None
    :type ph_min: float, optional
    :param ph_max: maximum value for the filter based on PH value, defaults to None
    :type ph_max: float, optional
    :param t_min: minimum value for the filter based on temperature value, defaults to None
    :type t_min: float, optional
    :param t_max: maximum value for the filter based on temperature value, defaults to None
    :type t_max: float, optional
    :return: tuple of lists (atom1_cs,atom2_cs as list)
    :rtype: tuple
    """

    x = []
    y = []
    cs_data = get_data_from_bmrb(residue=residue, ph_min=ph_min, ph_max=ph_max,
                                 t_min=t_min, t_max=t_max, standard_amino_acids=False)
    cs_dict = _list_to_dict(cs_data[0], cs_data[1])
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


def get_filtered_data_from_bmrb(residue: str,
                                atom: str,
                                filtering_rules: Optional[list],
                                ph_min: Optional[float] = None,
                                ph_max: Optional[float] = None,
                                t_min: Optional[float] = None,
                                t_max: Optional[float] = None,
                                h_tolerance: Optional[float] = 0.1,
                                c_tolerance: Optional[float] = 2.0,
                                n_tolerance: Optional[float] = 2.0,
                                standard_amino_acids: Optional[bool] = True) -> list:
    """
    Fetches the chemical shift data for a given residue and filters them using chemical shift values of other atoms in
    the same residue as a filter

    :param residue: residue name in IUPAC format
    :type residue: str
    :param atom: atom name in IUPAC format
    :type atom: str
    :param filtering_rules: list of atoms and chemical shift values as
        tuples to use as a filter; example[('CA',64.5),('H',7.8)]
    :type filtering_rules: list
    :param ph_min: minimum value for the filter based on PH value, defaults to None
    :type ph_min: float, optional
    :param ph_max: maximum value for the filter based on PH value, defaults to None
    :type ph_max: float, optional
    :param t_min: minimum value for the filter based on temperature value, defaults to None
    :type t_min: float, optional
    :param t_max: maximum value for the filter based on temperature value, defaults to None
    :type t_max: float, optional
    :param h_tolerance: tolerance value in ppm to filter H chemical shifts. It will allow the chemical shift value
        given in the filtering_rules +/- tolerance, defaults to 0.1
    :type h_tolerance: float, optional
    :param c_tolerance: tolerance value in ppm to filter C chemical shifts, It will allow the chemical shift value
        given in the filtering_rules +/- tolerance, defaults to 2.0
    :type c_tolerance: float, optional
    :param n_tolerance: tolerance value in ppm to filter N chemical shifts, It will allow the chemical shift value
        given in the filtering_rules +/- tolerance, defaults to 2.0
    :type n_tolerance: float, optional
    :param standard_amino_acids: get data only form standard amino acids and nucleic acids when '*' is used for residue,
        defaults to True
    :type standard_amino_acids: bool, optional
    :return: chemical shift values as list
    :rtype: list
    """

    cs_data = get_data_from_bmrb(residue=residue, ph_min=ph_min, ph_max=ph_max,
                                 t_min=t_min, t_max=t_max, standard_amino_acids=standard_amino_acids)

    def cs_filt(cs_dict2: dict, atm2: str, cs_val: float) -> dict:
        cs_width = 0.5
        if 'H' in atm2:
            cs_width = h_tolerance
        if 'C' in atm2:
            cs_width = c_tolerance
        if 'N' in atm2:
            cs_width = n_tolerance
        out_cs_dict = {}
        for key2 in cs_dict2.keys():
            try:
                if cs_val - cs_width <= cs_dict2[key2][atm2][0] <= cs_val + cs_width:
                    out_cs_dict[key2] = cs_dict2[key2]
            except KeyError:
                pass
        return out_cs_dict

    cs_dict = _list_to_dict(cs_data[0], cs_data[1])
    for rule in filtering_rules:
        cs_dict = cs_filt(cs_dict, rule[0], rule[1])
    x = []
    for key in cs_dict.keys():
        for atm in cs_dict[key].keys():
            if atm == atom:
                x.append(cs_dict[key][atm][0])
    return x


def get_statistics(residue: str,
                   atom: str,
                   list_of_atoms: Optional[Union[str, List[str]]] = None,
                   filtered: Optional[bool] = True,
                   sd_limit: Optional[float] = 10.0,
                   ambiguity: Optional[Union[str, int]] = '*',
                   ph_min: Optional[float] = None,
                   ph_max: Optional[float] = None,
                   t_min: Optional[float] = None,
                   t_max: Optional[float] = None,
                   standard_amino_acids: Optional[bool] = True,
                   verbose: bool = False) -> dict:
    """
    Provides chemical shift statistics like Mean, Median, Standard deviation for a given residue/atom/list of atoms

    :param residue: Residue name in IUPAC format; use '*' for any residue
    :type residue: str
    :param atom: Atom name in IUPAC format; Wild card supported; example HB*; use '*' for any atom
    :type atom: str
    :param list_of_atoms: list of atoms; example ['ALA-CB','CYS-N','TYR-CB'], defaults to None
    :type list_of_atoms: list, optional
    :param filtered: Filters values beyond (sd_limit)*(standard deviation) on both sides of the mean,
        defaults to True
    :type filtered: bool, optional
    :param sd_limit: scaling factor used to filter data based on standard deviation, defaults to 10.0
    :type sd_limit: float, optional
    :param ambiguity: filter based on chemical shift ambiguity code, defaults to '*' means no filter
    :type ambiguity: str, optional
    :param ph_min: minimum value for the filter based on PH value, defaults to None
    :type ph_min: float, optional
    :param ph_max: maximum value for the filter based on PH value, defaults to None
    :type ph_max: float, optional
    :param t_min: minimum value for the filter based on temperature value, defaults to None
    :type t_min: float, optional
    :param t_max: maximum value for the filter based on temperature value, defaults to None
    :type t_max: float, optional
    :param standard_amino_acids: get data only form standard amino acids and nucleic acids when '*' is used for residue,
        defaults to True
    :type standard_amino_acids: bool, optional
    :param verbose: display the statistics on screen; defaults  to False
    :type verbose: bool, optional
    :return: chemical shift statistics as a dictionary
    :rtype: dict
    """

    if residue is None and atom is None and list_of_atoms is None:
        logging.error('Please provide residue name or atom name or list of atoms as a list ')
        raise TypeError(
            'Please provide at least one of the three positional arguments: residue (or) atom (or) list_of_atoms')
    columns, data = get_data_from_bmrb(residue=residue,
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
        x = [float(i[cs_index]) for i in data if i[res_index] == res and i[atm_index] == at]
        cs_stat[atm]['cs'] = x
        cs_stat[atm]['mean'] = round(float(numpy.mean(x)), 3)
        cs_stat[atm]['median'] = round(float(numpy.median(x)), 3)
        cs_stat[atm]['std'] = round(float(numpy.std(x)), 3)
        cs_stat[atm]['min'] = min(x)
        cs_stat[atm]['max'] = max(x)
        cs_stat[atm]['count'] = len(x)
    if verbose:
        for key in cs_stat.keys():
            print('\t\t{}'.format(key))
            print('\t\tMean\t\t\t\t:{}'.format(cs_stat[key]['mean']))
            print('\t\tMedian\t\t\t\t:{}'.format(cs_stat[key]['median']))
            print('\t\tStandard deviation\t\t:{}'.format(cs_stat[key]['std']))
            print('\t\tMinimum\t\t\t\t:{}'.format(cs_stat[key]['min']))
            print('\t\tMaximum\t\t\t\t:{}'.format(cs_stat[key]['max']))
            print('\t\tCount\t\t\t\t:{}'.format(cs_stat[key]['count']))
            print("\n")
    return cs_stat

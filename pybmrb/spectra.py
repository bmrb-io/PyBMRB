#!/usr/bin/env python3


import logging
import csv
from pybmrb.get_entry_cs_data import ChemicalShift
import plotly.express as px


class Spectra(object):
    '''
    Converts one dimensional chemical shift list into multidimensional peak list and generates interactive plots
    '''

    @classmethod
    def create_c13hsqc_peaklist(self, bmrb_ids, input_file_names=None, auth_tag=False,
                                draw_trace=False):
        '''
        Converts one dimensional chemical shifts from list of BMRB entries into CHSQC peak list

        :param bmrb_ids: BMRB entry ID or list of entry ids
        :param input_file_names: Input NMR-STAR file name
        :param auth_tag: Use author provided sequence numbering from BMRB/NMR-STAR file; default: False
        :param draw_trace: Connect the matching residues using sequence numbering; default: False
        :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track);cs_track is a dictionary { matching atoms:[cs_values]}
        '''
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
        cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids, auth_tag=auth_tag)
        cs_data.update(cs_data_bmrb)
        if input_file_names is not None:
            cs_data_file = ChemicalShift.from_file(input_file_names, auth_tag=auth_tag)
            cs_data.update(cs_data_file)
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

    @classmethod
    def create_tocsy_peaklist(self, bmrb_ids, input_file_names=None, auth_tag=False,
                              draw_trace=False):
        """
        Converts one dimensional chemical shifts from list of BMRB entries into TOCSY peak list

        :param bmrb_ids: BMRB entry ID or list of entry ids
        :param input_file_names: Input NMR-STAR file name
        :param auth_tag: Use author provided sequence numbering from BMRB/NMR-STAR file; default: False
        :param draw_trace: Connect the matching residues using sequence numbering; default: False
        :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track);cs_track is a dictionary { matching atoms:[cs_values]}
        """
        cs_data = {}
        cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids, auth_tag=auth_tag)
        cs_data.update(cs_data_bmrb)
        if input_file_names is not None:
            cs_data_file = ChemicalShift.from_file(input_file_names, auth_tag=auth_tag)
            cs_data.update(cs_data_file)
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

    @classmethod
    def create_2d_peaklist(self, bmrb_ids, atom_x, atom_y, input_file_names=None, auth_tag=False,
                           draw_trace=False):
        """
         Converts one dimensional chemical shifts from list of BMRB entries into generic 2D peak list

        :param bmrb_ids: BMRB entry ID or list of entry ids
        :param atom_x: atom for x coordinate in IUPAC format
        :param atom_y: atom for y coordinate in IUPAC format
        :param input_file_names: Input NMR-STAR file name
        :param auth_tag: Use author provided sequence numbering from BMRB/NMR-STAR file; default: False
        :param draw_trace: Connect the matching residues using sequence numbering; default: False
        :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track);cs_track is a dictionary { matching atoms:[cs_values]}
        """
        cs_data = {}
        cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids, auth_tag=auth_tag)
        cs_data.update(cs_data_bmrb)
        if input_file_names is not None:
            cs_data_file = ChemicalShift.from_file(input_file_names, auth_tag=auth_tag)
            cs_data.update(cs_data_file)
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

    @classmethod
    def create_n15hsqc_peaklist(self, bmrb_ids, input_file_names=None, auth_tag=False, draw_trace=False,
                                include_sidechain=True):
        '''
        Converts one dimensional chemical shifts from list of BMRB entries into NHSQC peak list

        :param bmrb_ids: BMRB entry ID or list of entry ids
        :param input_file_names: Input NMR-STAR file name
        :param auth_tag: Use author provided sequence numbering from BMRB/NMR-STAR file; default: False
        :param draw_trace: Connect the matching residues using sequence numbering; default: False
        :param include_sidechain: include side chain NHs ; default: True
        :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track);cs_track is a dictionary { matching atoms:[cs_values]}
        '''
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
            cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids, auth_tag=auth_tag)
            cs_data.update(cs_data_bmrb)
        if input_file_names is not None:
            cs_data_file = ChemicalShift.from_file(input_file_names, auth_tag=auth_tag)
            cs_data.update(cs_data_file)
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
                        tag = '{}-{}-{}-{}'.format(data_id, chain, seq_no, residue)
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

    @classmethod
    def n15hsqc(self, bmrb_ids=None, input_file_names=None, auth_tag=False, legend=None, draw_trace=False,
                include_sidechain=True,
                peak_list=None,
                output_format='html',
                output_file=None,
                output_image_width=800,
                output_image_height=600,
                show_visualization=True):
        '''
        Plots NHSQC spectrum  for a given list of BMRB IDs (or) local NMR-STAR files (or) both

        :param bmrb_ids: list of BMRB IDs or single BMRB ID default None
        :param input_file_names: list of NMR-STAR files ; default None
        :param auth_tag: use author provided sequence numbering from BMRB /NMR-STAR file; default: False
        :param legend: legend based on residue name or data set id; values: None, residue, datase ; default None
        :param draw_trace: Connect matching residues by list
        :param include_sidechain: include side-chain NH peaks
        :param peak_list: optionally peak list as csv file cam be provided
        :param output_format: html,jpg,png,pdf,wabp,None; default None opens figure in default web browser
        :param output_file: output file name default None
        :param output_image_width: output image width; default 800
        :param output_image_height: output image height; default 600
        :param show_visualization: Automatically opens the visualization on a web browser; default True
        :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track);cs_track is a dictionary { matching atoms:[cs_values]}
        '''
        if peak_list is not None:
            x1 = []
            y1 = []
            with open(peak_list) as csvfile:
                spamreader = csv.reader(csvfile, delimiter=',')
                for row in spamreader:
                    x1.append(float(row[0]))
                    y1.append(float(row[1]))
        atom_x = 'H'
        atom_y = 'N'
        peak_list_2d = self.create_n15hsqc_peaklist(bmrb_ids,
                                                    input_file_names=input_file_names,
                                                    auth_tag=auth_tag,
                                                    draw_trace=draw_trace,
                                                    include_sidechain=include_sidechain)
        x = peak_list_2d[0]
        y = peak_list_2d[1]
        data_set = peak_list_2d[2]
        info = peak_list_2d[3]
        res = peak_list_2d[4]
        cs_track = peak_list_2d[5]
        if len(x) == 0 or len(y) == 0:
            logging.error('Resuired chemical shifts not found')
            raise ValueError('Required chemical shifts not found')
        if legend is None:
            fig = px.scatter(x=x, y=y,
                             symbol=data_set,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>15</sup>N (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x1, y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_layout(showlegend=False)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'residue':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     # "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>15</sup>N (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, mode='lines', opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'dataset':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=data_set,
                             labels={"color": "Data set",
                                     # "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>15</sup>N (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")
        if show_visualization: fig.show()
        if output_file is not None:
            if output_format == 'html':
                if output_file.split(".")[-1] == 'html':
                    fig.write_html('{}'.format(output_file))
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_html('{}.html'.format(output_file))
                    logging.info('Successfully written {}.html'.format(output_file))
            elif output_format == 'jpg':
                if output_file.split(".")[-1] == 'jpg':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.jpg'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.jpg'.format(output_file))
            elif output_format == 'png':
                if output_file.split(".")[-1] == 'png':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.png'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.png'.format(output_file))
            elif output_format == 'pdf':
                if output_file.split(".")[-1] == 'pdf':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.pdf'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.pdf'.format(output_file))
            elif output_format == 'webp':
                if output_file.split(".")[-1] == 'webp':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.webp'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.wepb'.format(output_file))
            else:
                logging.ERROR('Output file format nor support:{}'.format(output_format))
        return x, y, data_set, info, res, cs_track

    @classmethod
    def c13hsqc(self, bmrb_ids, input_file_names=None, auth_tag=False, legend=None, draw_trace=False,
                peak_list=None,
                output_format=None,
                output_file=None,
                output_image_width=800,
                output_image_height=600,
                show_visualization=True):
        '''
        Plots CHSQC spectrum  for a given list of BMRB IDs (or) local NMR-STAR files (or) both

        :param bmrb_ids: list of BMRB IDs or single BMRB ID default: None
        :param input_file_names: list of NMR-STAR files : default: None
        :param auth_tag: use author provided sequence numbering from BMRB /NMR-STAR file; default: False
        :param legend: legend based on residue name or data set id; values: None, residue, datase ; default None
        :param draw_trace: Connect matching residues by list
        :param peak_list: optionally peak list as csv file cam be provided
        :param output_format: html,jpg,png,pdf,wabp,None; default None opens figure in default web browser
        :param output_file: output file name default None
        :param output_image_width: output image width; default 800
        :param output_image_height: output image height; default 600
        :param show_visualization: Automatically opens the visualization on a web browser; default True
        :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track);cs_track is a dictionary { matching atoms:[cs_values]}
        '''
        peak_list_2d = self.create_c13hsqc_peaklist(bmrb_ids,
                                                    input_file_names=input_file_names,
                                                    auth_tag=auth_tag,
                                                    draw_trace=draw_trace)
        if peak_list is not None:
            x1 = []
            y1 = []
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
            logging.error('Resuired chemical shifts not found')
            raise ValueError('Required chemical shifts not found')
        if legend is None:
            fig = px.scatter(x=x, y=y,
                             symbol=data_set,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>13</sup>C (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_layout(showlegend=False)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'residue':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     # "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>13</sup>C (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, mode='lines', opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'dataset':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=data_set,
                             labels={"color": "Data set",
                                     # "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>13</sup>C (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")
        if show_visualization: fig.show()
        if output_file is not None:
            if output_format == 'html':
                if output_file.split(".")[-1] == 'html':
                    fig.write_html('{}'.format(output_file))
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_html('{}.html'.format(output_file))
                    logging.info('Successfully written {}.html'.format(output_file))
            elif output_format == 'jpg':
                if output_file.split(".")[-1] == 'jpg':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.jpg'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.jpg'.format(output_file))
            elif output_format == 'png':
                if output_file.split(".")[-1] == 'png':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.png'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.png'.format(output_file))
            elif output_format == 'pdf':
                if output_file.split(".")[-1] == 'pdf':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.pdf'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.pdf'.format(output_file))
            elif output_format == 'webp':
                if output_file.split(".")[-1] == 'webp':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.webp'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.wepb'.format(output_file))
            else:
                logging.ERROR('Output file format nor support:{}'.format(output_format))
        return x, y, data_set, info, res, cs_track

    @classmethod
    def tocsy(self, bmrb_ids, input_file_names=None, auth_tag=False, legend=None, draw_trace=False,
              peak_list=None,
              output_format=None,
              output_file=None,
              output_image_width=800,
              output_image_height=600,
              show_visualization=True):
        '''
        Plots TOCSY spectrum  for a given list of BMRB IDs (or) local NMR-STAR files (or) both

        :param bmrb_ids: list of BMRB IDs or single BMRB ID default: None
        :param input_file_names: list of NMR-STAR files : default: None
        :param auth_tag: use author provided sequence numbering from BMRB /NMR-STAR file; default: False
        :param legend: legend based on residue name or data set id; values: None, residue, datase ; default None
        :param draw_trace: Connect matching residues by list
        :param peak_list: optionally peak list as csv file cam be provided
        :param output_format: html,jpg,png,pdf,wabp,None; default None opens figure in default web browser
        :param output_file: output file name default None
        :param output_image_width: output image width; default 800
        :param output_image_height: output image height; default 600
        :param show_visualization: Automatically opens the visualization on a web browser; default True
        :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track);cs_track is a dictionary { matching atoms:[cs_values]}
        '''
        peak_list_2d = self.create_tocsy_peaklist(bmrb_ids,
                                                  input_file_names=input_file_names,
                                                  auth_tag=auth_tag,
                                                  draw_trace=draw_trace)
        if peak_list is not None:
            x1 = []
            y1 = []
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
            logging.error('Resuired chemical shifts not found')
            raise ValueError('Required chemical shifts not found')
        if legend is None:
            fig = px.scatter(x=x, y=y,
                             symbol=data_set,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>1</sup>H (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_layout(showlegend=False)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'residue':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     # "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>1</sup>H (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, mode='lines', opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'dataset':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=data_set,
                             labels={"color": "Data set",
                                     # "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>1</sup>H (ppm)'}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")
        if show_visualization: fig.show()
        if output_file is not None:
            if output_format == 'html':
                if output_file.split(".")[-1] == 'html':
                    fig.write_html('{}'.format(output_file))
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_html('{}.html'.format(output_file))
                    logging.info('Successfully written {}.html'.format(output_file))
            elif output_format == 'jpg':
                if output_file.split(".")[-1] == 'jpg':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.jpg'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.jpg'.format(output_file))
            elif output_format == 'png':
                if output_file.split(".")[-1] == 'png':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.png'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.png'.format(output_file))
            elif output_format == 'pdf':
                if output_file.split(".")[-1] == 'pdf':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.pdf'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.pdf'.format(output_file))
            elif output_format == 'webp':
                if output_file.split(".")[-1] == 'webp':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.webp'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.wepb'.format(output_file))
            else:
                logging.ERROR('Output file format nor support:{}'.format(output_format))
        return x, y, data_set, info, res, cs_track

    @classmethod
    def generic_2d(self, bmrb_ids, input_file_names=None, atom_x='H', atom_y='N', auth_tag=False, legend=None,
                   draw_trace=False,
                   peak_list=None,
                   output_format=None,
                   output_file=None,
                   output_image_width=800,
                   output_image_height=600,
                   show_visualization=True):
        '''
        Plots generic 2D spectrum  for a given list of BMRB IDs (or) local NMR-STAR files (or) both

        :param bmrb_ids: list of BMRB IDs or single BMRB ID default: None
        :param input_file_names: list of NMR-STAR files : default: None
        :param atom_x: atom name for x coordinate in IUPAC format; default : 'H
        :param atom_y: atom name for y coordinate in IUPAC format; default :'N'
        :param auth_tag: use author provided sequence numbering from BMRB /NMR-STAR file; default: False
        :param legend: legend based on residue name or data set id; values: None, residue, datase ; default None
        :param draw_trace: Connect matching residues by list
        :param include_sidechain: include side-chain NH peaks
        :param peak_list: optionally peak list as csv file cam be provided
        :param output_format: html,jpg,png,pdf,wabp,None; default None opens figure in default web browser
        :param output_file: output file name default None
        :param output_image_width: output image width; default 800
        :param output_image_height: output image height; default 600
        :param show_visualization: Automatically opens the visualization on a web browser; default True
        :return: tuple of lists and dictionary (x,y,data_set,info,res,cs_track);cs_track is a dictionary { matching atoms:[cs_values]}
        '''
        peak_list_2d = self.create_2d_peaklist(bmrb_ids, atom_x=atom_x, atom_y=atom_y,
                                               input_file_names=input_file_names,
                                               auth_tag=auth_tag,
                                               draw_trace=draw_trace)
        if peak_list is not None:
            x1 = []
            y1 = []
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
            logging.error('Resuired chemical shifts not found')
            raise ValueError('Required chemical shifts not found')
        if legend is None:
            fig = px.scatter(x=x, y=y,
                             symbol=data_set,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     "symbol": "Data set",
                                     "x": '{} (ppm)'.format(atom_x),
                                     "y": '{} (ppm)'.format(atom_y)}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_layout(showlegend=False)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'residue':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     # "symbol": "Data set",
                                     "x": '{} (ppm)'.format(atom_x),
                                     "y": '{} (ppm)'.format(atom_y)}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, mode='lines', opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'dataset':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=data_set,
                             labels={"color": "Data set",
                                     # "symbol": "Data set",
                                     "x": '{} (ppm)'.format(atom_x),
                                     "y": '{} (ppm)'.format(atom_y)}, opacity=0.7)
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, opacity=0.7)
            if peak_list is not None:
                fig.add_scatter(x=x1, y=y1, mode='markers', name='Peak list', opacity=0.7)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")
        if show_visualization: fig.show()
        if output_file is not None:
            if output_format == 'html':
                if output_file.split(".")[-1] == 'html':
                    fig.write_html('{}'.format(output_file))
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_html('{}.html'.format(output_file))
                    logging.info('Successfully written {}.html'.format(output_file))
            elif output_format == 'jpg':
                if output_file.split(".")[-1] == 'jpg':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.jpg'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.jpg'.format(output_file))
            elif output_format == 'png':
                if output_file.split(".")[-1] == 'png':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.png'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.png'.format(output_file))
            elif output_format == 'pdf':
                if output_file.split(".")[-1] == 'pdf':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.pdf'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.pdf'.format(output_file))
            elif output_format == 'webp':
                if output_file.split(".")[-1] == 'webp':
                    fig.write_image('{}'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}'.format(output_file))
                else:
                    fig.write_image('{}.webp'.format(output_file), width=output_image_width, height=output_image_height)
                    logging.info('Successfully written {}.wepb'.format(output_file))
            else:
                logging.ERROR('Output file format nor support:{}'.format(output_format))
        return x, y, data_set, info, res, cs_track

# if __name__ == "__main__":
# # Generating examples for documentation
# Spectra.n15hsqc(bmrb_ids=15060,output_format='jpg',legend='residue',output_file='../docs/_images/15060_n15',
#                 show_visualization=False)
# Spectra.n15hsqc(bmrb_ids=15060, output_format='html', legend='residue',output_file='../docs/_static/15060_n15',
#                 show_visualization=False)
# Spectra.n15hsqc(bmrb_ids=[17076,17077], file_names='tests/test_data/MyData.str',
#                 output_format='jpg', output_file='../docs/_images/multi_n15',
#                 legend='dataset',show_visualization=False)
# Spectra.n15hsqc(bmrb_ids=[17076,17077], file_names='tests/test_data/MyData.str',
#                 output_format='html', output_file='../docs/_static/multi_n15',
#                 legend='dataset',show_visualization=False)
# Spectra.n15hsqc(bmrb_ids=[ 17076, 17077], file_names='tests/test_data/MyData.str',
#                 peak_list='tests/test_data/test_peak_list.csv',
#                 output_format='jpg', output_file='../docs/_images/multi_n152',
#                 legend='dataset', show_visualization=False)
# Spectra.n15hsqc(bmrb_ids=[ 17076, 17077], file_names='tests/test_data/MyData.str',
#                 peak_list='tests/test_data/test_peak_list.csv',
#                 output_format='html', output_file='../docs/_static/multi_n152',
#                 legend='dataset', show_visualization=False)
# Spectra.n15hsqc(bmrb_ids=[ 17076, 17077],file_names='tests/test_data/MyData.str',
#                 output_format='jpg', output_file='../docs/_images/multi2_n15',
#                 legend='dataset',draw_trace=True, show_visualization=False,)
# Spectra.n15hsqc(bmrb_ids=[ 17076, 17077], file_names='tests/test_data/MyData.str',
#                 output_format='html', output_file='../docs/_static/multi2_n15',
#                 legend='dataset',draw_trace=True,show_visualization=False)
#
# Spectra.c13hsqc(bmrb_ids=15060, output_format='jpg', legend='residue', output_file='../docs/_images/15060_c13',
#                 show_visualization=False)
# Spectra.c13hsqc(bmrb_ids=15060, output_format='html', legend='residue', output_file='../docs/_static/15060_c13',
#                 show_visualization=False)
# Spectra.c13hsqc(bmrb_ids=[17074, 17076, 17077], output_format='jpg', output_file='../docs/_images/multi_c13',
#                 legend='dataset', show_visualization=False)
# Spectra.c13hsqc(bmrb_ids=[17074, 17076, 17077], output_format='html', output_file='../docs/_static/multi_c13',
#                 legend='dataset', show_visualization=False)
# Spectra.n15hsqc(bmrb_ids=[17074, 17076, 17077], output_format='jpg', output_file='../docs/_images/multi2_c13',
#                 legend='dataset', draw_trace=True, show_visualization=False, )
# Spectra.c13hsqc(bmrb_ids=[17074, 17076, 17077], output_format='html', output_file='../docs/_static/multi2_c13',
#                 legend='dataset', draw_trace=True, show_visualization=False)
#
# Spectra.tocsy(bmrb_ids=15060, output_format='jpg', legend='residue', output_file='../docs/_images/15060_tocsy',
#                 show_visualization=False)
# Spectra.tocsy(bmrb_ids=15060, output_format='html', legend='residue', output_file='../docs/_static/15060_tocsy',
#                 show_visualization=False)
# Spectra.tocsy(bmrb_ids=[17074, 17076, 17077], output_format='jpg', output_file='../docs/_images/multi_tocsy',
#                 legend='dataset', show_visualization=False)
# Spectra.tocsy(bmrb_ids=[17074, 17076, 17077], output_format='html', output_file='../docs/_static/multi_tocsy',
#                 legend='dataset', show_visualization=False)
# Spectra.tocsy(bmrb_ids=[17074, 17076, 17077], output_format='jpg', output_file='../docs/_images/multi_tocsy2',
#               legend='residue', show_visualization=False)
# Spectra.tocsy(bmrb_ids=[17074, 17076, 17077], output_format='html', output_file='../docs/_static/multi_tocsy2',
#               legend='residue', show_visualization=False)
# Spectra.tocsy(bmrb_ids=[17074, 17076, 17077], output_format='jpg', output_file='../docs/_images/multi2_tocsy',
#                 legend='dataset', draw_trace=True, show_visualization=False, )
# Spectra.tocsy(bmrb_ids=[17074, 17076, 17077], output_format='html', output_file='../docs/_static/multi2_tocsy',
#                 legend='dataset', draw_trace=True, show_visualization=False)
#
# Spectra.generic_2d(bmrb_ids=15060, atom_x='N',atom_y='CB', output_format='jpg', legend='residue',
#                    output_file='../docs/_images/15060_2d',
#                    show_visualization=False)
# Spectra.generic_2d(bmrb_ids=15060, atom_x='N',atom_y='CB',output_format='html', legend='residue',
#                    output_file='../docs/_static/15060_2d', show_visualization=False)
# Spectra.generic_2d(bmrb_ids=[17074, 17076, 17077], atom_x='N',atom_y='CB',output_format='jpg',
#                    output_file='../docs/_images/multi_2d', legend='dataset', show_visualization=False)
# Spectra.generic_2d(bmrb_ids=[17074, 17076, 17077], atom_x='N',atom_y='CB',output_format='html',
#                    output_file='../docs/_static/multi_2d', legend='dataset', show_visualization=False)
# Spectra.generic_2d(bmrb_ids=[17074, 17076, 17077], atom_x='N',atom_y='CB',output_format='jpg',
#                    output_file='../docs/_images/multi2_2d', legend='dataset', draw_trace=True,
#                    show_visualization=False, )
# Spectra.generic_2d(bmrb_ids=[17074, 17076, 17077], atom_x='N',atom_y='CB',output_format='html',
#                    output_file='../docs/_static/multi2_2d',
#                    legend='dataset', draw_trace=True, show_visualization=False)

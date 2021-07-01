import logging

from pybmrb2.get_cs_data import ChemicalShift
import plotly.express as px
import plotly.graph_objects as go

class Spectra(object):

    def __init__(self):
        pass

    @classmethod
    def create_c13hsqc_peaklist(self, bmrb_ids, file_names=None, atom_x='H', atom_y='N', auth_tag=False,
                                draw_trace=False):
        cs_data = {}
        cs_data_bmrb = ChemicalShift.from_bmrb(bmrb_ids, auth_tag=auth_tag)
        cs_data.update(cs_data_bmrb)
        if file_names is not None:
            cs_data_file = ChemicalShift.from_file(file_names, auth_tag=auth_tag)
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
                    except KeyError:
                        logging.debug('Data not found:{},{},{}'.format(data_id, chain, seq_no))
            cs_track = {}
            if draw_trace:
                for k in atom_ids.keys():
                    if len(atom_ids[k][0]) > 1:
                        cs_track[k] = atom_ids[k]
        return x, y, data_set, info, res, cs_track

    @classmethod
    def create_n15hsqc_peaklist(self, bmrb_ids, file_names=None, atom_x='H', atom_y='N', auth_tag=False, draw_trace=False,
                                include_sidechain=False):
        sidechain_nh_atoms = {'ARG': {
            'ARG-HH11-NH1': ['HH11','NH1'],
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
                'HIS-HD1-ND1': ['HD1','ND1'],
                'HIS-HE2-ND1': ['HE2', 'NE2']},
            'TRP': {
                'TRP-HE1-NE1': ['HE1', 'NE1']},
            'LYS': {
                'LYS-HZ-NZ': ['HZ', 'NZ'],
                'LYS-HZ1-NZ': ['HZ1', 'NZ'],
                'LYS-HZ2-NZ': ['HZ2', 'NZ'],
                'LYS-HZ3-NZ': ['HZ3', 'NZ']}
        }
        cs_data={}
        cs_data_bmrb=ChemicalShift.from_bmrb(bmrb_ids,auth_tag=auth_tag)
        cs_data.update(cs_data_bmrb)
        if file_names is not None:
            cs_data_file=ChemicalShift.from_file(file_names,auth_tag=auth_tag)
            cs_data.update(cs_data_file)
        data_set=[]
        x=[]
        y=[]
        res=[]
        info=[]
        atom_ids={}
        for data_id in cs_data.keys():
            for chain in cs_data[data_id].keys():
                for seq_no in cs_data[data_id][chain]['seq_ids']:
                    try:
                        x_cs=cs_data[data_id][chain][seq_no][atom_x][2]
                        y_cs=cs_data[data_id][chain][seq_no][atom_y][2]
                        residue=cs_data[data_id][chain][seq_no][atom_y][0]
                        res.append(residue)
                        tag='{}-{}-{}-{}'.format(data_id,chain,seq_no,residue)
                        data_set.append(data_id)
                        atom_id='{}-{}-{}'.format(chain,seq_no,residue)
                        if draw_trace:
                            if atom_id not in atom_ids.keys():
                                atom_ids[atom_id]=[[],[]]
                            atom_ids[atom_id][0].append(x_cs)
                            atom_ids[atom_id][1].append(y_cs)
                        x.append(x_cs)
                        y.append(y_cs)
                        info.append(tag)
                        if include_sidechain:
                            if residue in sidechain_nh_atoms.keys():
                                for atom_list in sidechain_nh_atoms[residue]:
                                    ax=sidechain_nh_atoms[residue][atom_list][0]
                                    ay=sidechain_nh_atoms[residue][atom_list][1]
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
                        logging.debug('Data not found:{},{},{}'.format(data_id,chain,seq_no))
            cs_track={}
            if draw_trace:
                for k in atom_ids.keys():
                    if len(atom_ids[k][0])>1:
                        cs_track[k]=atom_ids[k]
        return x,y,data_set,info,res,cs_track

    @classmethod
    def n15hsqc(self,bmrb_ids,file_names=None,auth_tag=False,legend = None,draw_trace=False,
                include_sidechain=True,
                output_format='html',
                output_file=None,
                output_image_width=800,
                output_image_height=600):
        atom_x='H'
        atom_y='N'
        peak_list_2d = self.create_n15hsqc_peaklist(bmrb_ids,
                                                    file_names=file_names,
                                                    atom_x=atom_x,
                                                    atom_y=atom_y,
                                                    auth_tag=auth_tag,
                                                    draw_trace=draw_trace,
                                                    include_sidechain=include_sidechain)
        x = peak_list_2d[0]
        y = peak_list_2d[1]
        data_set = peak_list_2d[2]
        info = peak_list_2d[3]
        res = peak_list_2d[4]
        cs_track = peak_list_2d[5]
        if legend is None:
            fig = px.scatter(x=x, y=y,
                             symbol=data_set,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>15</sup>N (ppm)'})
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k)
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
                                     "y": '<sup>15</sup>N (ppm)'})
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k, mode='lines')
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend == 'dataset':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=data_set,
                             labels={"color": "Data set",
                                     # "symbol": "Data set",
                                     "x": '<sup>1</sup>H (ppm)',
                                     "y": '<sup>15</sup>N (ppm)'})
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")
        fig.show()
        if output_file is not None:
            if output_format=='html':
                fig.write_html('{}.html'.format(output_file))
                logging.info('Sucessfully written {}.html'.format(output_file))
            elif output_format=='jpg':
                fig.write_image('{}.jpg'.format(output_file),width=output_image_width,height=output_image_height)
                logging.info('Sucessfully written {}.jpg'.format(output_file))
            elif output_format=='png':
                fig.write_image('{}.png'.format(output_file),width=output_image_width,height=output_image_height)
                logging.info('Sucessfully written {}.png'.format(output_file))
            elif output_format=='pdf':
                fig.write_image('{}.pdf'.format(output_file),width=output_image_width,height=output_image_height)
                logging.info('Sucessfully written {}.pdf'.format(output_file))
            elif output_format=='webp':
                fig.write_image('{}.webp'.format(output_file),width=output_image_width,height=output_image_height)
                logging.info('Sucessfully written {}.wepb'.format(output_file))
            else:
                logging.ERROR('Output file format nor support:{}'.format(output_format))
        return True

    def generic_2d(self,bmrb_ids,file_names=None,atom_x='H',atom_y='N',auth_tag=False,legend = None,draw_trace=False):

        peak_list_2d=self.create_n15hsqc_peaklist(bmrb_ids,
                                                  file_names=file_names,
                                                  atom_x=atom_x,
                                                  atom_y=atom_y,
                                                  auth_tag=auth_tag,
                                                  draw_trace=draw_trace)
        x=peak_list_2d[0]
        y=peak_list_2d[1]
        data_set=peak_list_2d[2]
        info=peak_list_2d[3]
        res=peak_list_2d[4]
        cs_track=peak_list_2d[5]

        if legend is None:
            fig=px.scatter(x=x,y=y,
                           symbol=data_set,
                           hover_name=info,
                           color=res,
                           labels={"color":"Residue",
                                   "symbol":"Data set",
                                   "x":'{} (ppm)'.format(atom_x),
                                   "y":'{} (ppm)'.format(atom_y)})
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0], y=cs_track[k][1], name=k)
            fig.update_layout(showlegend=False)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend=='residue':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=res,
                             labels={"color": "Residue",
                                     #"symbol": "Data set",
                                     "x": '{} (ppm)'.format(atom_x),
                                     "y": '{} (ppm)'.format(atom_y)})
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0],y=cs_track[k][1],name=k,mode='lines')
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")

        elif legend=='dataset':
            fig = px.scatter(x=x, y=y,
                             hover_name=info,
                             color=data_set,
                             labels={"color": "Data set",
                                     #"symbol": "Data set",
                                     "x": '{} (ppm)'.format(atom_x),
                                     "y": '{} (ppm)'.format(atom_y)})
            if draw_trace:
                for k in cs_track.keys():
                    fig.add_scatter(x=cs_track[k][0],y=cs_track[k][1],name=k)
            fig.update_xaxes(autorange="reversed")
            fig.update_yaxes(autorange="reversed")
        fig.show()







if __name__ == "__main__":
    Spectra.n15hsqc([17074,17077,17076],legend='residue',output_file='test',output_format='webp')
    #s.generic_2d([15060,18857],file_names='/Users/Kumaran/MyData.str',draw_trace=True,legend='dataset')

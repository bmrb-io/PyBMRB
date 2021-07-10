#!/usr/bin/env python3

import logging
import plotly.express as px
from pybmrb.get_database_cs_data import ChemicalShiftStatistics

class Histogram(object):
    '''
    Fetches chemical shift records from BMRB, plots histograms and calculates statistical properties.
    '''


    @classmethod
    def hist(self, residue=None, atom=None, list_of_atoms=None, filtered=True, sd_limit=10, ambiguity='*',
                           ph_min=None, ph_max=None, t_min=None, t_max=None,
             histnorm='',
             standard_amino_acids=True,
             plot_type='histogram', output_format='html',
             output_file=None,
             output_image_width=800,
             output_image_height=600,
             show_visualization = True
        ):
        '''
        plots histogram for a given list of atoms and residues with some filters. One of either residue or atom or list of atoms is required

        :param residue:  residue name in IUPAC format; '*' for all standard residues: default None
        :param atom:  atom name in IUPAC format; '*' for all standard atoms; default None
        :param list_of_atoms: list of atoms in IUPAC actom; example '['ALA-CA','CYS-CB']'; default None
        :param filtered: Filters values beyond (sd_limt)*(standard deviation) on both sides of the mean; default:True
        :param sd_limit: scaling factor used to filter data based on standard deviation; default 10
        :param ambiguity: ambiguity filter; default '*' => no filter
        :param ph_min: PH filter (min);default None
        :param ph_max: PH filter (max); default None
        :param t_min: Temperature filter (min); default None
        :param histnrom: histnorm for the distribution 'probability','percent','probability density')
        :param t_max: Temperature filter (max); default None
        :param standard_amino_acids: get data only form 20 natural amino acids,4 standard DNA and 4 standard RNA; default:True
        :param plot_type: plot type; supported types 'histogram','box','violin' ; default histogram
        :param output_format: output format type; supported types 'html','jpg','png','pdf','webp';default 'html'
        :param output_file: output file name; if provided, the output will be written in a file , otherwise opens is a web browser; default ;None
        :param output_image_width: output image width to write in a file; default:800
        :param output_image_height: output image height to write in a file; default 600
        :param show_visualization: Automatically opens the visualization on a web browser; default True
        :return: chemical shift data and tags as tuple (chemical shifts, tags)
        '''
        columns, cs_data = ChemicalShiftStatistics.get_data_from_bmrb(residue=residue,
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
        cs_index = columns.index('Atom_chem_shift.Val')
        ph_index = columns.index('Sample_conditions.pH')
        temp_index = columns.index('Sample_conditions.Temperature_K')
        res_index = columns.index('Atom_chem_shift.Comp_ID')
        atm_index = columns.index('Atom_chem_shift.Atom_ID')
        amb_index = columns.index('Atom_chem_shift.Ambiguity_code')
        if len(cs_data)==0:
            logging.error('No matching atom, or values found in the database')
            raise ValueError('No matching atom, or values found in the database')
        x=[]
        tag=[]
        for row in cs_data:
            t='{}-{}'.format(row[res_index],row[atm_index])
            tag.append(t)
            x.append(row[cs_index])
        if len(x)==0:
            logging.error('No matching atom, or values found in the database')
            raise ValueError('No matching atom, or values found in the database')
        if plot_type=='histogram':
            fig = px.histogram(x,color=tag,histnorm=histnorm,
                               labels={"color": "Atom",
                                         "value": 'Chemical shift (ppm)',
                                         "count": 'Count'},opacity=0.5)
            fig.update_xaxes(autorange="reversed")
            fig.update_layout(barmode='overlay')
        elif plot_type=='box':
            fig = px.box(x=tag,y=x,color=tag,
                               labels={"color": "Atom",
                                       "x":"",
                                       "y":'Chemical shift (ppm)'
                                         })
            fig.update_xaxes(tickangle=90)
        elif plot_type=='violin':
            fig=px.violin(x=tag,y=x,color=tag,
                               labels={"color": "Atom",
                                         "x":"",
                                       "y":'Chemical shift (ppm)'
                                         })
            fig.update_xaxes(tickangle=90)
        else:
            logging.error('Plot type not supported : {}'.format(plot_type))
            raise TypeError('Plot type not supported : {}'.format(plot_type))
        if show_visualization: fig.show()
        if output_file is not None:
            if output_format == 'html':
                fig.write_html('{}.html'.format(output_file))
                logging.info('Sucessfully written {}.html'.format(output_file))
            elif output_format == 'jpg':
                fig.write_image('{}.jpg'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.jpg'.format(output_file))
            elif output_format == 'png':
                fig.write_image('{}.png'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.png'.format(output_file))
            elif output_format == 'pdf':
                fig.write_image('{}.pdf'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.pdf'.format(output_file))
            elif output_format == 'webp':
                fig.write_image('{}.webp'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.wepb'.format(output_file))
            else:
                logging.ERROR('Output file format nor support:{}'.format(output_format))
        return x,tag


    @classmethod
    def hist2d(cls,residue,atom1,atom2,filtered=True, sd_limit=10,
                               ambiguity1='*', ambiguity2='*',
                               ph_min=None,ph_max=None,t_min=None,t_max=None,
               histnorm='',
               plot_type='heatmap',output_format='html',
                output_file=None,
                output_image_width=800,
                output_image_height=600,
               show_visualization=True):
        '''
        Generates chemical shift correlation plot for any two atoms from a given residue.

        :param residue: residue name in IUPAC format
        :param atom1: atom name in IUPAC format
        :param atom2: atom name in IUPAC format
        :param filtered: Filters values beyond (sd_limt)*(standard deviation) on both sides of the mean; default:True
        :param sd_limit: scaling factor used to filter data based on standard deviation; default 10
        :param ambiguity: ambiguity filter; default '*' => no filter
        :param ph_min: PH filter (min);default None
        :param ph_max: PH filter (max); default None
        :param t_min: Temperature filter (min); default None
        :param t_max: Temperature filter (max); default None
        :param histnrom: histnorm for the distribution 'probability','percent','probability density')
        :param plot_type: plot type; support types 'heatmap','contour'
        :param output_format: output format type; supported types 'html','jpg','png','pdf','webp';default 'html'
        :param output_file: output file name; if provided, the output will be written in a file , otherwise opens is a web browser; default ;None
        :param output_image_width: output image width to write in a file; default:800
        :param output_image_height: output image height to write in a file; default 600
        :param show_visualization: Automatically opens the visualization on a web browser; default True
        :return: tuple (chemical shift list of atom1, chemical shift list of atom2)
        '''
        x,y = ChemicalShiftStatistics.get_2d_chemical_shifts(residue=residue,
                                                             atom1=atom1,
                                                             atom2=atom2,
                                                             filtered=filtered,
                                                             sd_limit=sd_limit,
                                                             ambiguity1=ambiguity1,
                                                             ambiguity2=ambiguity2,
                                                             ph_min=ph_min,
                                                             ph_max=ph_max,
                                                             t_min=t_min,
                                                             t_max=t_max)
        if len(x)==0:
            logging.error('No matching atom, or values found in the database')
            raise ValueError('No matching atom, or values found in the database')
        if plot_type == 'heatmap':
            fig = px.density_heatmap(x=x,y=y, marginal_x="histogram", marginal_y="histogram",
                                     histnorm=histnorm,
                                     labels={
                                             "x": '{} (ppm)'.format(atom1),
                                             "y": '{} (ppm)'.format(atom2)},
                                     )
            fig.update_layout(xaxis=dict(autorange='reversed'),
                              yaxis=dict(autorange='reversed'),
                              xaxis2=dict(showticklabels=True),
                              yaxis3=dict(showticklabels=True))
        elif plot_type == 'contour':
            fig = px.density_contour(x=x, y=y, marginal_x="histogram", marginal_y="histogram",histnorm=histnorm,
                                     labels={
                                         "x": '{} (ppm)'.format(atom1),
                                         "y": '{} (ppm)'.format(atom2)},
                                     )
            fig.update_layout(xaxis=dict(autorange='reversed'),
                              yaxis=dict(autorange='reversed'),
                              xaxis2=dict(showticklabels=True),
                              yaxis3=dict(showticklabels=True))
        else:
            logging.error('Plot type not supported : {}'.format(plot_type))
            raise TypeError('Plot type not supported : {}'.format(plot_type))
        if show_visualization: fig.show()
        if output_file is not None:
            if output_format == 'html':
                fig.write_html('{}.html'.format(output_file))
                logging.info('Sucessfully written {}.html'.format(output_file))
            elif output_format == 'jpg':
                fig.write_image('{}.jpg'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.jpg'.format(output_file))
            elif output_format == 'png':
                fig.write_image('{}.png'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.png'.format(output_file))
            elif output_format == 'pdf':
                fig.write_image('{}.pdf'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.pdf'.format(output_file))
            elif output_format == 'webp':
                fig.write_image('{}.webp'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.wepb'.format(output_file))
            else:
                logging.ERROR('Output file format nor support:{}'.format(output_format))
        return x,y

    @classmethod
    def conditional_hist(cls, residue, atom, filtering_rules,
                         ph_min=None, ph_max=None, t_min=None, t_max=None,
                         histnorm='',
                         standard_amino_acids=True,
                         plot_type='histogram', output_format='html',
                         output_file=None,
                         output_image_width=800,
                         output_image_height=600,
                         show_visualization=True
                         ):
        '''
        Plots the distribution of the given atom in the residue along with the filtered distribution besed
        on the chemical shift values of the other atoms in the residue

        :param residue: residue name in IUPAC format; example 'CYS'
        :param atom: atom name in IUPAC format; example 'CB'
        :param filtering_rules: list of atoms and chemical shift values as tuples; example[('CA',64.5),('H',7.8)]
        :param ph_min: PH filter (min);default None
        :param ph_max: PH filter (max); default None
        :param t_min: Temperature filter (min); default None
        :param t_max: Temperature filter (max); default None
        :param histnrom: histnorm for the distribution 'probability','percent','probability density')
        :param standard_amino_acids: get data only form 20 natural amino acids,4 standard DNA and 4 standard RNA; default:True
        :param plot_type: plot type; support types 'heatmap','contour'
        :param output_format: output format type; supported types 'html','jpg','png','pdf','webp';default 'html'
        :param output_file: output file name; if provided, the output will be written in a file , otherwise opens is a web browser; default ;None
        :param output_image_width: output image width to write in a file; default:800
        :param output_image_height: output image height to write in a file; default 600
        :param show_visualization: Automatically opens the visualization on a web browser; default True
        :return: chemical shift data and tags as tuple (chemical shifts, tags)
        '''
        columns, cs_data = ChemicalShiftStatistics.get_data_from_bmrb(residue=residue,
                                                                      atom=atom,
                                                                      ph_min=ph_min,
                                                                      ph_max=ph_max,
                                                                      t_min=t_min,
                                                                      t_max=t_max,
                                                                      standard_amino_acids=standard_amino_acids)
        cs_index = columns.index('Atom_chem_shift.Val')
        ph_index = columns.index('Sample_conditions.pH')
        temp_index = columns.index('Sample_conditions.Temperature_K')
        res_index = columns.index('Atom_chem_shift.Comp_ID')
        atm_index = columns.index('Atom_chem_shift.Atom_ID')
        amb_index = columns.index('Atom_chem_shift.Ambiguity_code')
        if len(cs_data)==0:
            logging.error('No matching atom, or values found in the database')
            raise ValueError('No matching atom, or values found in the database')
        x1=ChemicalShiftStatistics.get_filtered_data_from_bmrb(residue=residue,
                                                              atom=atom,
                                                              filtering_rules=filtering_rules,
                                                              ph_min=ph_min,
                                                              ph_max=ph_max,
                                                              t_min=t_min,
                                                              t_max=t_max,
                                                              standard_amino_acids=standard_amino_acids
                                                              )
        if len(x1)==0:
            logging.error('No matching atom, or values found in the database')
            raise ValueError('No matching atom, or values found in the database')

        x = []
        tag = []
        for row in cs_data:
            t = '{}-{}'.format(row[res_index], row[atm_index])
            tag.append(t)
            x.append(row[cs_index])
        for i in x1:
            x.append(i)
            tag.append('Filtered')
        if len(x)==0:
            logging.error('No matching atom, or values found in the database')
            raise ValueError('No matching atom, or values found in the database')
        if plot_type == 'histogram':
            fig = px.histogram(x, color=tag,histnorm=histnorm,
                               labels={"color": "Atom",
                                       "value": 'Chemical shift (ppm)',
                                       "count": 'Count'}, opacity=0.5)
            fig.update_xaxes(autorange="reversed")
        elif plot_type == 'box':
            fig = px.box(x=tag, y=x, color=tag,
                         labels={"color": "Atom",
                                 "x": "",
                                 "y": 'Chemical shift (ppm)'
                                 })
            fig.update_xaxes(tickangle=90)
        elif plot_type == 'violin':
            fig = px.violin(x=tag, y=x, color=tag,
                            labels={"color": "Atom",
                                    "x": "",
                                    "y": 'Chemical shift (ppm)'
                                    })
            fig.update_xaxes(tickangle=90)
        else:
            logging.error('Plot type not supported : {}'.format(plot_type))
            raise TypeError('Plot type not supported : {}'.format(plot_type))
        fig.update_layout(barmode='overlay')
        if show_visualization: fig.show()
        if output_file is not None:
            if output_format == 'html':
                fig.write_html('{}.html'.format(output_file))
                logging.info('Sucessfully written {}.html'.format(output_file))
            elif output_format == 'jpg':
                fig.write_image('{}.jpg'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.jpg'.format(output_file))
            elif output_format == 'png':
                fig.write_image('{}.png'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.png'.format(output_file))
            elif output_format == 'pdf':
                fig.write_image('{}.pdf'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.pdf'.format(output_file))
            elif output_format == 'webp':
                fig.write_image('{}.webp'.format(output_file), width=output_image_width, height=output_image_height)
                logging.info('Sucessfully written {}.wepb'.format(output_file))
            else:
                logging.ERROR('Output file format nor support:{}'.format(output_format))
        return x, tag

if __name__=="__main__":
    Histogram.hist(list_of_atoms=['GLN-CB','CYS-CB','TYR-CB'],histnorm='probability density',
                   output_format='jpg',output_file='../docs/_images/multi_hist')
    Histogram.hist(list_of_atoms=['GLN-CB', 'CYS-CB', 'TYR-CB'], histnorm='probability density',
                   output_format='html', output_file='../docs/_static/multi_hist')
    Histogram.hist(list_of_atoms=['GLN-CB', 'CYS-CB', 'TYR-CB'], plot_type='violin',
                   output_format='jpg', output_file='../docs/_images/multi_violin')
    Histogram.hist(list_of_atoms=['GLN-CB', 'CYS-CB', 'TYR-CB'], plot_type='violin',
                   output_format='html', output_file='../docs/_static/multi_violin')

    # Histogram.hist(residue='CYS',atom='CB',output_format='jpg',
    #                output_file='../docs/_images/cys_cb_hist',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='html',
    #                output_file='../docs/_static/cys_cb_hist',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='jpg',
    #                sd_limit=5, output_file='../docs/_images/cys_cb_hist_sd5',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='html',
    #                sd_limit=5, output_file='../docs/_static/cys_cb_hist_sd5',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='jpg',
    #                ph_min=7.0,ph_max=8.2, output_file='../docs/_images/cys_cb_hist_ph',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='html',
    #                ph_min=7.0,ph_max=8.2, output_file='../docs/_static/cys_cb_hist_ph',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='jpg', plot_type='box',
    #                 output_file='../docs/_images/cys_cb_box_sd5',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='html', plot_type='box',
    #                output_file='../docs/_static/cys_cb_box_sd5',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='jpg', plot_type='violin',
    #                 output_file='../docs/_images/cys_cb_violin_sd5',show_visualization=False)
    # Histogram.hist(residue='CYS', atom='CB', output_format='html', plot_type='violin',
    #                 output_file='../docs/_static/cys_cb_violin_sd5',show_visualization=False)
    # Histogram.hist(residue='GLN', atom='H*', output_format='jpg',
    #                histnorm='probability density',
    #                output_file='../docs/_images/gln_h_hist',show_visualization=False)
    # Histogram.hist(residue='GLN', atom='H*', output_format='html',
    #                histnorm='probability density',
    #                 output_file='../docs/_static/gln_h_hist',show_visualization=False)
    # Histogram.hist(residue='ASP', atom='*', output_format='jpg',
    #                output_file='../docs/_images/asp_hist',show_visualization=False)
    # Histogram.hist(residue='ASP', atom='*', output_format='html',
    #                output_file='../docs/_static/asp_hist',show_visualization=False)
    # Histogram.hist(residue='*', atom='CG*', output_format='jpg',histnorm='percent',
    #                output_file='../docs/_images/cg_hist',show_visualization=False)
    # Histogram.hist(residue='*', atom='CG*', output_format='html', histnorm='percent',
    #                output_file='../docs/_static/cg_hist',show_visualization=False)
    #
    #
    # Histogram.hist2d(residue='CYS',atom1='CA',atom2='CB',output_format='jpg', sd_limit=5,
    #                  output_file='../docs/_images/cys-ca-cb',show_visualization=False)
    # Histogram.hist2d(residue='CYS', atom1='CA', atom2='CB', output_format='html', sd_limit=5,
    #                  output_file='../docs/_static/cys-ca-cb',show_visualization=False)
    # Histogram.hist2d(residue='GLN',atom1='HE21',atom2='HE22',output_format='jpg', sd_limit=5,
    #                  plot_type='contour', output_file='../docs/_images/gln-2d',show_visualization=False)
    # Histogram.hist2d(residue='GLN', atom1='HE21', atom2='HE22', output_format='html', sd_limit=5,
    #                  plot_type='contour',output_file='../docs/_static/gln-2d',show_visualization=False)
    #
    # Histogram.conditional_hist(residue='CYS',atom='CB', histnorm='percent',
    #                            filtering_rules=[('H',8.9)],output_format='jpg',
    #                            output_file='../docs/_images/filt1',show_visualization=False)
    # Histogram.conditional_hist(residue='CYS', atom='CB', histnorm='percent',
    #                            filtering_rules=[('H', 8.9)], output_format='html',
    #                            output_file='../docs/_static/filt1',show_visualization=False)
    # Histogram.conditional_hist(residue='CYS', atom='CB', histnorm='percent',
    #                            filtering_rules=[('H', 8.9), ('CA', 61)], output_format='jpg',
    #                            output_file='../docs/_images/filt2',show_visualization=False)
    # Histogram.conditional_hist(residue='CYS', atom='CB', histnorm='percent',
    #                            filtering_rules=[('H', 8.9), ('CA', 61)],output_format='html',
    #                            output_file='../docs/_static/filt2',show_visualization=False)
    #
    #
    #

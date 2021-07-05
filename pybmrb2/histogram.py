import logging
import csv
import plotly.express as px
from pybmrb2.get_database_cs_data import ChemicalShiftStatistics

class Histogram(object):

    def __init__(self):
        pass
    @classmethod
    def hist(self, residue=None, atom=None, list_of_atoms=None, filtered=True, sd_limit=10, ambiguity='*',
                           ph_min=None, ph_max=None, t_min=None, t_max=None,
             standard_amino_acids=True,
             plot_type='histogram', output_format='html',
             output_file=None,
             output_image_width=800,
             output_image_height=600
        ):
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

        x=[]
        tag=[]
        for row in cs_data:
            t='{}-{}'.format(row[res_index],row[atm_index])
            tag.append(t)
            x.append(row[cs_index])
        if plot_type=='histogram':
            fig = px.histogram(x,color=tag,
                               labels={"color": "Atom",
                                         "value": 'Chemical shift (ppm)',
                                         "count": 'Count'},opacity=0.5)
            fig.update_xaxes(autorange="reversed")
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
        fig.show()
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


    @classmethod
    def hist2d(cls,residue,atom1,atom2,filtered=True, sd_limit=10,
                               ambiguity1='*', ambiguity2='*',
                               ph_min=None,ph_max=None,t_min=None,t_max=None,
               plot_type='heatmap',output_format='html',
                output_file=None,
                output_image_width=800,
                output_image_height=600):
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
        if plot_type == 'heatmap':
            fig = px.density_heatmap(x=x,y=y, marginal_x="histogram", marginal_y="histogram",
                                     labels={
                                             "x": '{} (ppm)'.format(atom1),
                                             "y": '{} (ppm)'.format(atom2)},
                                     )
            fig.update_layout(xaxis=dict(autorange='reversed'),
                              yaxis=dict(autorange='reversed'),
                              xaxis2=dict(showticklabels=True),
                              yaxis3=dict(showticklabels=True))
        elif plot_type == 'contour':
            fig = px.density_contour(x=x, y=y, marginal_x="histogram", marginal_y="histogram",
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
        fig.show()
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

if __name__=="__main__":
    Histogram.hist(residue='CYS',atom='H',sd_limit=10)




import logging
import csv
import plotly.express as px
from pybmrb2.get_database_cs_data import ChemicalShiftStatistics

class Histogram(object):

    def __init__(self):
        pass
    @classmethod
    def hist(self, residue=None, atom=None, list_of_atoms=None, filtered=True, sd_limit=10, ambiguity='*',
                           ph_min=None, ph_max=None, t_min=None, t_max=None, standard_amino_acids=True):
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

        fig = px.histogram(x,color=tag,
                           labels={"color": "Atom",
                                     "value": 'Chemical shift (ppm)',
                                     "count": 'Count'},opacity=0.5)
        fig.update_xaxes(autorange="reversed")
        fig.show()

    @classmethod
    def hist2d(cls,residue,atom1,atom2,filtered=True, sd_limit=10,
                               ambiguity1='*', ambiguity2='*',
                               ph_min=None,ph_max=None,t_min=None,t_max=None,standard_amino_acids=True):
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
        fig = px.density_heatmap(x=x,y=y, marginal_x="histogram", marginal_y="histogram",
                                 labels={
                                         "x": '{} (ppm)'.format(atom1),
                                         "y": '{} (ppm)'.format(atom2)},
                                 )
        fig.update_xaxes(autorange="reversed")
        fig.update_yaxes(autorange="reversed")
        fig.show()

if __name__=="__main__":
    Histogram.hist(list_of_atoms=['ALA-CB','TYR-N','GLY-N'])




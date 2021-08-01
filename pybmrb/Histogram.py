#!/usr/bin/env python3
"""
This module generates visualizations for chemical shift distribution of a given list of atoms and residues. It can also
generate chemical shift correlation plots between two atoms in the same residue and filtered chemical shift
distribution based on user defined rules
"""
import logging
import plotly.express as px
from pybmrb import ChemicalShiftStatistics
from typing import Union, List, Optional


def hist(residue: Optional[str] = None,
         atom: Optional[str] = None,
         list_of_atoms: Optional[Union[str, List[str]]] = None,
         filtered: Optional[bool] = True,
         sd_limit: Optional[float] = 10.0,
         ambiguity: Optional[Union[str, int]] = '*',
         ph_min: Optional[float] = None,
         ph_max: Optional[float] = None,
         t_min: Optional[float] = None,
         t_max: Optional[float] = None,
         standard_amino_acids: Optional[bool] = True,
         histnorm: Optional[str] = None,
         plot_type: Optional[str] = 'histogram',
         output_format: Optional[str] = 'html',
         output_file: Optional[str] = None,
         output_image_width: Optional[int] = 800,
         output_image_height: Optional[int] = 600,
         show_visualization: Optional[bool] = True
         ) -> tuple:
    """
    plots histogram for a given list of atoms and residues with some filters. One of either residue or atom or list
    of atoms is required

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
    :param histnorm: Specifies the type of normalization used for this histogram trace. If None, the span of each bar
        corresponds to the number of occurrences (i.e. the number of data points lying inside the bins).
        If 'percent' / 'probability', the span of each bar corresponds to the percentage / fraction of occurrences with
        respect to the total number of sample points (here, the sum of all bin HEIGHTS equals 100% / 1). If 'density',
        the span of each bar corresponds to the number of occurrences in a bin divided by the size of the bin interval
        (here, the sum of all bin AREAS equals the total number of sample points). If 'probability density', the area of
        each bar corresponds to the probability that an event will fall into the corresponding bin (here, the sum of all
        bin AREAS equals 1) defaults to None
    :type histnorm: str, optional
    :param plot_type: visualization type 'histogram', 'box' or 'violin', defaults to 'histogram'
    :type plot_type: str, optional
    :param output_format: visualizations can be exported as interactive 'html' file
        or as static images in 'jpg','jpeg','png','pdf','webp','svg', defaults to 'html'
    :type output_format: str, optional
    :param output_file: file name to export visualization
    :type output_file: str, optional
    :param output_image_width: The width of the exported image in layout pixels, defaults to 800
    :type output_image_width: int, optional
    :param output_image_height: The height of the exported image in layout pixels, defaults to 600
    :type output_image_height: int, optional
    :param show_visualization: Visualization automatically opens in a web browser or as
        embedded visualization in Jupyter Notebook. This feature can be disabled
        by setting this flag as False, defaults to True
    :type show_visualization: bool, optional
    :return: chemical shift data and tags as tuple (chemical shifts, tags)
    :rtype: tuple
    """

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
    # ph_index = columns.index('Sample_conditions.pH')
    # temp_index = columns.index('Sample_conditions.Temperature_K')
    res_index = columns.index('Atom_chem_shift.Comp_ID')
    atm_index = columns.index('Atom_chem_shift.Atom_ID')
    # amb_index = columns.index('Atom_chem_shift.Ambiguity_code')
    if len(cs_data) == 0:
        logging.error('No matching atom, or values found in the database')
        raise ValueError('No matching atom, or values found in the database')
    x = []
    tag = []
    for row in cs_data:
        t = '{}-{}'.format(row[res_index], row[atm_index])
        tag.append(t)
        x.append(row[cs_index])
    if len(x) == 0:
        logging.error('No matching atom, or values found in the database')
        raise ValueError('No matching atom, or values found in the database')
    if plot_type == 'histogram':
        fig = px.histogram(x, color=tag, histnorm=histnorm,
                           title='Chemical shift statistics from BMRB',
                           labels={"color": "Atom",
                                   "value": 'Chemical shift (ppm)'
                                   }, opacity=0.5).update_layout(yaxis_title='Count')
        fig.update_xaxes(autorange="reversed")
        fig.update_layout(barmode='overlay')
        fig.update(layout=dict(title=dict(x=0.5)))
    elif plot_type == 'box':
        fig = px.box(x=tag, y=x, color=tag,
                     title='Chemical shift statistics from BMRB',
                     labels={"color": "Atom",
                             "x": "",
                             "y": 'Chemical shift (ppm)'
                             }).update(layout=dict(title=dict(x=0.5)))
        fig.update_xaxes(tickangle=90)
    elif plot_type == 'violin':
        fig = px.violin(x=tag, y=x, color=tag,
                        title='Chemical shift statistics from BMRB',
                        labels={"color": "Atom",
                                "x": "",
                                "y": 'Chemical shift (ppm)'
                                }).update(layout=dict(title=dict(x=0.5)))
        fig.update_xaxes(tickangle=90)
    else:
        logging.error('Plot type not supported : {}'.format(plot_type))
        raise TypeError('Plot type not supported : {}'.format(plot_type))
    if show_visualization:
        fig.show()
    if output_file is not None:
        if output_format == 'html':
            fig.write_html(output_file)
        elif output_format in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'json']:
            fig.write_image(file=output_file, format=output_format, width=output_image_width,
                            height=output_image_height)
        else:
            logging.error('Output file format not supported:{}'.format(output_format))
    return x, tag


def hist2d(residue: str,
           atom1: str,
           atom2: str,
           filtered: Optional[bool] = True,
           sd_limit: Optional[float] = 10.0,
           ambiguity1: Optional[Union[str, int]] = '*',
           ambiguity2: Optional[Union[str, int]] = '*',
           ph_min: Optional[float] = None,
           ph_max: Optional[float] = None,
           t_min: Optional[float] = None,
           t_max: Optional[float] = None,
           histnorm: Optional[str] = None,
           plot_type: Optional[str] = 'heatmap',
           output_format: Optional[str] = 'html',
           output_file: Optional[str] = None,
           output_image_width: Optional[int] = 800,
           output_image_height: Optional[int] = 600,
           show_visualization: Optional[bool] = True) -> tuple:
    """
    Generates chemical shift correlation plot for any two atoms from a given residue.

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
    :param histnorm: Specifies the type of normalization used for this histogram trace. If None, the span of each bar
        corresponds to the number of occurrences (i.e. the number of data points lying inside the bins).
        If 'percent' / 'probability', the span of each bar corresponds to the percentage / fraction of occurrences with
        respect to the total number of sample points (here, the sum of all bin HEIGHTS equals 100% / 1). If 'density',
        the span of each bar corresponds to the number of occurrences in a bin divided by the size of the bin interval
        (here, the sum of all bin AREAS equals the total number of sample points). If 'probability density', the area of
        each bar corresponds to the probability that an event will fall into the corresponding bin (here, the sum of all
        bin AREAS equals 1) defaults to None
    :type histnorm: str, optional
    :param plot_type: visualization type 'contour' or 'heatmap', defaults to 'heatmap'
    :type plot_type: str, optional
    :param output_format: visualizations can be exported as interactive 'html' file
        or as static images in 'jpg','jpeg','png','pdf','webp','svg', defaults to 'html'
    :type output_format: str, optional
    :param output_file: file name to export visualization
    :type output_file: str, optional
    :param output_image_width: The width of the exported image in layout pixels, defaults to 800
    :type output_image_width: int, optional
    :param output_image_height: The height of the exported image in layout pixels, defaults to 600
    :type output_image_height: int, optional
    :param show_visualization: Visualization automatically opens in a web browser or as
        embedded visualization in Jupyter Notebook. This feature can be disabled
        by setting this flag as False, defaults to True
    :type show_visualization: bool, optional
    :return: tuple (chemical shift list of atom1, chemical shift list of atom2)
    :rtype: tuple
    """

    x, y = ChemicalShiftStatistics.get_2d_chemical_shifts(residue=residue,
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
    if len(x) == 0:
        logging.error('No matching atom, or values found in the database')
        raise ValueError('No matching atom, or values found in the database')
    if plot_type == 'heatmap':
        fig = px.density_heatmap(x=x, y=y, marginal_x="histogram", marginal_y="histogram",
                                 title='Chemical shift statistics from BMRB',
                                 histnorm=histnorm,
                                 labels={
                                     "x": '{} (ppm)'.format(atom1),
                                     "y": '{} (ppm)'.format(atom2)},
                                 ).update(layout=dict(title=dict(x=0.5)))
        fig.update_layout(xaxis=dict(autorange='reversed'),
                          yaxis=dict(autorange='reversed'),
                          xaxis2=dict(showticklabels=True),
                          yaxis3=dict(showticklabels=True))
    elif plot_type == 'contour':
        fig = px.density_contour(x=x, y=y, marginal_x="histogram", marginal_y="histogram", histnorm=histnorm,
                                 title='Chemical shift statistics from BMRB',
                                 labels={
                                     "x": '{} (ppm)'.format(atom1),
                                     "y": '{} (ppm)'.format(atom2)},
                                 ).update(layout=dict(title=dict(x=0.5)))
        fig.update_layout(xaxis=dict(autorange='reversed'),
                          yaxis=dict(autorange='reversed'),
                          xaxis2=dict(showticklabels=True),
                          yaxis3=dict(showticklabels=True))
    else:
        logging.error('Plot type not supported : {}'.format(plot_type))
        raise TypeError('Plot type not supported : {}'.format(plot_type))
    if show_visualization:
        fig.show()
    if output_file is not None:
        if output_format == 'html':
            fig.write_html(output_file)
        elif output_format in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'json']:
            fig.write_image(file=output_file, format=output_format, width=output_image_width,
                            height=output_image_height)
        else:
            logging.error('Output file format not supported:{}'.format(output_format))
    return x, y


def conditional_hist(residue: str,
                     atom: str,
                     filtering_rules: Optional[list],
                     ph_min: Optional[float] = None,
                     ph_max: Optional[float] = None,
                     t_min: Optional[float] = None,
                     t_max: Optional[float] = None,
                     h_tolerance: Optional[float] = 0.1,
                     c_tolerance: Optional[float] = 2.0,
                     n_tolerance: Optional[float] = 2.0,
                     standard_amino_acids: Optional[bool] = True,
                     histnorm: Optional[str] = None,
                     plot_type: Optional[str] = 'histogram',
                     output_format: Optional[str] = 'html',
                     output_file: Optional[str] = None,
                     output_image_width: Optional[int] = 800,
                     output_image_height: Optional[int] = 600,
                     show_visualization: Optional[bool] = True) -> tuple:
    """
    Plots the distribution of the given atom in the residue along with the filtered distribution based
    on the chemical shift values of the other atoms in the residue

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
    :param histnorm: Specifies the type of normalization used for this histogram trace. If None, the span of each bar
        corresponds to the number of occurrences (i.e. the number of data points lying inside the bins).
        If 'percent' / 'probability', the span of each bar corresponds to the percentage / fraction of occurrences with
        respect to the total number of sample points (here, the sum of all bin HEIGHTS equals 100% / 1). If 'density',
        the span of each bar corresponds to the number of occurrences in a bin divided by the size of the bin interval
        (here, the sum of all bin AREAS equals the total number of sample points). If 'probability density', the area of
        each bar corresponds to the probability that an event will fall into the corresponding bin (here, the sum of all
        bin AREAS equals 1) defaults to None
    :type histnorm: str, optional
    :param plot_type: visualization type 'histogram', 'box' or 'violin', defaults to 'histogram'
    :type plot_type: str, optional
    :param output_format: visualizations can be exported as interactive 'html' file
        or as static images in 'jpg','jpeg','png','pdf','webp','svg', defaults to 'html'
    :type output_format: str, optional
    :param output_file: file name to export visualization
    :type output_file: str, optional
    :param output_image_width: The width of the exported image in layout pixels, defaults to 800
    :type output_image_width: int, optional
    :param output_image_height: The height of the exported image in layout pixels, defaults to 600
    :type output_image_height: int, optional
    :param show_visualization: Visualization automatically opens in a web browser or as
        embedded visualization in Jupyter Notebook. This feature can be disabled
        by setting this flag as False, defaults to True
    :type show_visualization: bool, optional
    :return: chemical shift data and tags as tuple (chemical shifts, tags)
    :rtype: tuple
    """
    columns, cs_data = ChemicalShiftStatistics.get_data_from_bmrb(residue=residue,
                                                                  atom=atom,
                                                                  ph_min=ph_min,
                                                                  ph_max=ph_max,
                                                                  t_min=t_min,
                                                                  t_max=t_max,
                                                                  standard_amino_acids=standard_amino_acids)
    cs_index = columns.index('Atom_chem_shift.Val')
    # ph_index = columns.index('Sample_conditions.pH')
    # temp_index = columns.index('Sample_conditions.Temperature_K')
    res_index = columns.index('Atom_chem_shift.Comp_ID')
    atm_index = columns.index('Atom_chem_shift.Atom_ID')
    # amb_index = columns.index('Atom_chem_shift.Ambiguity_code')
    if len(cs_data) == 0:
        logging.error('No matching atom, or values found in the database')
        raise ValueError('No matching atom, or values found in the database')
    x1 = ChemicalShiftStatistics.get_filtered_data_from_bmrb(residue=residue,
                                                             atom=atom,
                                                             filtering_rules=filtering_rules,
                                                             ph_min=ph_min,
                                                             ph_max=ph_max,
                                                             t_min=t_min,
                                                             t_max=t_max,
                                                             h_tolerance=h_tolerance,
                                                             c_tolerance=c_tolerance,
                                                             n_tolerance=n_tolerance,
                                                             standard_amino_acids=standard_amino_acids
                                                             )
    if len(x1) == 0:
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
    if len(x) == 0:
        logging.error('No matching atom, or values found in the database')
        raise ValueError('No matching atom, or values found in the database')
    if plot_type == 'histogram':
        fig = px.histogram(x, color=tag, histnorm=histnorm,
                           title='Chemical shift statistics from BMRB',
                           labels={"color": "Atom",
                                   "value": 'Chemical shift (ppm)',
                                   "count": 'Count'}, opacity=0.5).update(layout=dict(title=dict(x=0.5)))
        fig.update_xaxes(autorange="reversed")
    elif plot_type == 'box':
        fig = px.box(x=tag, y=x, color=tag,
                     title='Chemical shift statistics from BMRB',
                     labels={"color": "Atom",
                             "x": "",
                             "y": 'Chemical shift (ppm)'
                             }).update(layout=dict(title=dict(x=0.5)))
        fig.update_xaxes(tickangle=90)
    elif plot_type == 'violin':
        fig = px.violin(x=tag, y=x, color=tag,
                        title='Chemical shift statistics from BMRB',
                        labels={"color": "Atom",
                                "x": "",
                                "y": 'Chemical shift (ppm)'
                                }).update(layout=dict(title=dict(x=0.5)))
        fig.update_xaxes(tickangle=90)
    else:
        logging.error('Plot type not supported : {}'.format(plot_type))
        raise TypeError('Plot type not supported : {}'.format(plot_type))
    fig.update_layout(barmode='overlay')
    if show_visualization:
        fig.show()
    if output_file is not None:
        if output_format == 'html':
            fig.write_html(output_file)
        elif output_format in ['png', 'jpg', 'jpeg', 'webp', 'svg', 'pdf', 'json']:
            fig.write_image(file=output_file, format=output_format, width=output_image_width,
                            height=output_image_height)
        else:
            logging.error('Output file format not supported:{}'.format(output_format))
    return x, tag

# if __name__=="__main__":


# Following script generates example images for readthedocs
# hist(residue='TYR',atom='CB',
#      output_format='jpg',
#      output_file='../docs/_images/quick_star_hist1.jpg',
#      show_visualization=False)
# hist(residue='TYR', atom='CB',
#      output_format='html',
#      output_file='../docs/_static/quick_star_hist1.html',
#      show_visualization=False)
#
# hist(residue='CYS', atom='CB',
#      output_format='jpg',
#      plot_type='box',
#      output_file='../docs/_images/quick_star_hist2.jpg',
#      show_visualization=False)
# hist(residue='CYS', atom='CB',
#      output_format='html',
#      plot_type='box',
#      output_file='../docs/_static/quick_star_hist2.html',
#      show_visualization=False)
#
# hist(residue='CYS', atom='CB',
#      output_format='jpg',
#      plot_type='violin',
#      output_file='../docs/_images/quick_star_hist3.jpg',
#      show_visualization=False)
# hist(residue='CYS', atom='CB',
#      output_format='html',
#      plot_type='violin',
#      output_file='../docs/_static/quick_star_hist3.html',
#      show_visualization=False)
#
# hist(residue='TYR', atom='H*',
#      output_format='jpg',
#      output_file='../docs/_images/quick_star_hist4.jpg',
#      show_visualization=False)
# hist(residue='TYR', atom='H*',
#      output_format='html',
#      output_file='../docs/_static/quick_star_hist4.html',
#      show_visualization=False)
#
# hist( atom='CB',
#      output_format='jpg',
#      output_file='../docs/_images/quick_star_hist5.jpg',
#      show_visualization=False)
# hist( atom='CB',
#      output_format='html',
#      output_file='../docs/_static/quick_star_hist5.html',
#      show_visualization=False)
#
# hist2d(residue='CYS',
#        atom1='N',
#        atom2='CB',
#        sd_limit=5,
#        output_format='jpg',
#        output_file='../docs/_images/quick_star_hist6.jpg',
#        show_visualization=False)
# hist2d(residue='CYS',
#        atom1='N',
#        atom2='CB',
#        sd_limit=5,
#        output_format='html',
#        output_file='../docs/_static/quick_star_hist6.html',
#        show_visualization=False)
# hist(atom='CB*',output_format='jpg',output_file='../docs/_images/sample_cbhist.jpg',
#      show_visualization=False)
# hist(atom='CB*', output_format='html', output_file='../docs/_static/sample_cbhist.html',
#      show_visualization=False)
#
# hist2d(residue='CYS',atom1='N',atom2='CB',sd_limit=5,
#        output_file='../docs/_images/sample_hist2d.jpg',output_format='jpg',
#        show_visualization=False
#        )
# hist2d(residue='CYS', atom1='N', atom2='CB', sd_limit=5,output_format='html',
#        output_file='../docs/_static/sample_hist2d.html',
#        show_visualization=False
#        )
#
# hist(residue='CYS',
#      atom='CB',
#      output_format='jpg',
#      output_file='../docs/_images/example20.jpg',
#      show_visualization=False)
# hist(residue='CYS',
#      atom='CB',
#      output_format='html',
#      output_file='../docs/_static/example20.html',
#      show_visualization=False)
#
# hist(residue='CYS',
#      atom='CB',
#      sd_limit=5,
#      output_format='jpg',
#      output_file='../docs/_images/example21.jpg',
#      show_visualization=False)
# hist(residue='CYS',
#      atom='CB',
#      sd_limit=5,
#      output_format='html',
#      output_file='../docs/_static/example21.html',
#      show_visualization=False)
#
# hist(residue='CYS',
#      atom='CB',
#      sd_limit=5,
#      ph_min=7.0,
#      ph_max=8.2,
#      output_format='jpg',
#      output_file='../docs/_images/example22.jpg',
#      show_visualization=False)
# hist(residue='CYS',
#      atom='CB',
#      sd_limit=5,
#      ph_min=7.0,
#      ph_max=8.2,
#      output_format='html',
#      output_file='../docs/_static/example22.html',
#      show_visualization=False)
#
# hist(residue='CYS',
#      atom='CB',
#      plot_type='box',
#      output_format='jpg',
#      output_file='../docs/_images/example23.jpg',
#      show_visualization=False)
# hist(residue='CYS',
#      atom='CB',
#      plot_type='box',
#      output_format='html',
#      output_file='../docs/_static/example23.html',
#      show_visualization=False)
#
# hist(residue='CYS',
#      atom='CB',
#      plot_type='violin',
#      output_format='jpg',
#      output_file='../docs/_images/example24.jpg',
#      show_visualization=False)
# hist(residue='CYS',
#      atom='CB',
#      plot_type='violin',
#      output_format='html',
#      output_file='../docs/_static/example24.html',
#      show_visualization=False)
#
# hist(list_of_atoms=['GLN-CB','CYS-CB','TYR-CB'],
#      histnorm='probability density',
#      output_format='jpg',
#      output_file='../docs/_images/example25.jpg',
#      show_visualization=False)
# hist(list_of_atoms=['GLN-CB', 'CYS-CB', 'TYR-CB'],
#      histnorm='probability density',
#      output_format='html',
#      output_file='../docs/_static/example25.html',
#      show_visualization=False)
#
# hist(list_of_atoms=['GLN-CB','CYS-CB','TYR-CB'],
#      plot_type='violin',
#      output_format='jpg',
#      output_file='../docs/_images/example26.jpg',
#      show_visualization=False)
# hist(list_of_atoms=['GLN-CB', 'CYS-CB', 'TYR-CB'],
#      plot_type='violin',
#      output_format='html',
#      output_file='../docs/_static/example26.html',
#      show_visualization=False)
#
# hist(residue='GLN', atom='H*', histnorm='probability density',
#      output_format='jpg',
#      output_file='../docs/_images/example27.jpg',
#      show_visualization=False)
# hist(residue='GLN', atom='H*', histnorm='probability density',
#      output_format='html',
#      output_file='../docs/_static/example27.html',
#      show_visualization=False)
#
# hist(residue='ASP',
#      sd_limit=5,
#      output_format='jpg',
#      output_file='../docs/_images/example28.jpg',
#      show_visualization=False)
# hist(residue='ASP',
#      sd_limit=5.0,
#      output_format='html',
#      output_file='../docs/_static/example28.html',
#      show_visualization=False)
#
# hist(atom='CG*',
#      histnorm='percent',
#      output_format='jpg',
#      output_file='../docs/_images/example29.jpg',
#      show_visualization=False)
# hist(atom='CG*',
#      histnorm='percent',
#      output_format='html',
#      output_file='../docs/_static/example29.html',
#      show_visualization=False)
#
# hist2d(residue='CYS', atom1='CA', atom2='CB', sd_limit=5,
#        output_format='jpg',
#        output_file='../docs/_images/example30.jpg',
#        show_visualization=False
#        )
# hist2d(residue='CYS', atom1='CA', atom2='CB', sd_limit=5,
#        output_format='html',
#        output_file='../docs/_static/example30.html',
#        show_visualization=False
#        )
#
# hist2d(residue='GLN', atom1='HE21', atom2='HE22', sd_limit=5,
#        output_format='jpg',
#        plot_type='contour',
#        output_file='../docs/_images/example31.jpg',
#        show_visualization=False
#        )
# hist2d(residue='GLN', atom1='HE21', atom2='HE22', sd_limit=5,
#        output_format='html',
#        plot_type='contour',
#        output_file='../docs/_static/example31.html',
#        show_visualization=False
#        )

# conditional_hist(residue='CYS',
#                  atom='CB', histnorm='percent',
#                  filtering_rules=[('H', 8.9)],
#                  output_format='jpg',
#                  output_file='../docs/_images/example32.jpg',
#                  show_visualization=False
#                  )
# conditional_hist(residue='CYS',
#                  atom='CB', histnorm='percent',
#                  filtering_rules=[('H', 8.9)],
#                  output_format='html',
#                  output_file='../docs/_static/example32.html',
#                  show_visualization=False
#                  )
#
# conditional_hist(residue='CYS',
#                  atom='CB', histnorm='percent',
#                  filtering_rules=[('H', 8.9),('CA', 61)],
#                  output_format='jpg',
#                  output_file='../docs/_images/example33.jpg',
#                  show_visualization=False
#                  )
# conditional_hist(residue='CYS',
#                  atom='CB', histnorm='percent',
#                  filtering_rules=[('H', 8.9),('CA', 61)],
#                  output_format='html',
#                  output_file='../docs/_static/example33.html',
#                  show_visualization=False
#                  )

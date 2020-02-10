from __future__ import print_function

import json
import ntpath
import optparse
import os
import sys
import re

import numpy as np
import plotly
import pynmrstar

# Determine if we are running in python3
PY3 = (sys.version_info[0] == 3)
(scriptPath, scriptName) = ntpath.split(os.path.realpath(__file__))
# pylint: disable=wrong-import-position,no-name-in-module
# pylint: disable=import-error,wrong-import-order
# Python version dependent loads
if PY3:
    from urllib.request import urlopen, Request
else:
    from urllib2 import urlopen, Request

NOTEBOOK = False
_OPACITY = 1.0
_AUTOOPEN = True
__version__ = "1.2.8"


class LACS_Viz(object):

    def __init__(self):
        d = self.read_lacs_output('data/2m07_1_A.lacs_out')
        self.plot_data(d)

    def read_lacs_output(self, file_name):
        dat = open(file_name, 'r').read()
        data = re.findall('\s+(\d+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\d+)\n', dat)
        off_atoms = [i for i in re.findall('OFFATOMS:([\S \t]+)\n', dat)[0].split(" ") if i != ""]
        off_values = [i for i in re.findall('OFFVALUES:\s*([\S \t]+)\n', dat)[0].split(" ") if i != ""]
        print(off_atoms, off_values)
        data_sets = []
        off_atoms = ['HN' if i == 'H' else i for i in off_atoms]
        off_atoms = ['NH' if i == 'N' else i for i in off_atoms]
        off_atoms = ['CO' if i == 'C' else i for i in off_atoms]
        for atm in off_atoms:
            data_dict = {}
            x_normal = []
            y_normal = []
            txt_normal = []
            x_outlier = []
            y_outlier = []
            txt_outlier = []
            x_line1 = []
            y_line1 = []
            txt_line1 = []
            x_line2 = []
            y_line2 = []
            txt_line2 = []
            for i in data:
                if i[3] == atm:
                    if int(i[6]) == 1:
                        x_normal.append(float(i[4]))
                        y_normal.append(float(i[5]))
                        txt_normal.append('{}-{}'.format(i[0], i[1]))
                    if int(i[6]) == 0:
                        x_outlier.append(float(i[4]))
                        y_outlier.append(float(i[5]))
                        txt_outlier.append('{}-{}'.format(i[0], i[1]))
                    if int(i[6]) == 2:
                        x_line1.append(float(i[4]))
                        y_line1.append(float(i[5]))
                        txt_line1.append('line 1')
                    if int(i[6]) == 3:
                        x_line2.append(float(i[4]))
                        y_line2.append(float(i[5]))
                        txt_line2.append('line 2')
            if len(x_normal):
                data_dict['normal'] = (x_normal, y_normal, txt_normal)
                data_dict['outlier'] = (x_outlier, y_outlier, txt_outlier)
                data_dict['line1'] = (x_line1, y_line1, txt_line1)
                data_dict['line2'] = (x_line2, y_line2, txt_line1)
                data_dict['atom'] = atm
                data_dict['offset'] = off_values[off_atoms.index(atm)]
                data_sets.append(data_dict)
        return data_sets

    def plot_data(self, data_sets):

        for dat in data_sets:
            data = []
            data.append(plotly.graph_objs.Scatter(x=dat['normal'][0],
                                                  y=dat['normal'][1],
                                                  text=dat['normal'][2],
                                                  mode='markers',
                                                  opacity=_OPACITY,
                                                  name='Normal'))
            data.append(plotly.graph_objs.Scatter(x=dat['outlier'][0],
                                                  y=dat['outlier'][1],
                                                  text=dat['outlier'][2],
                                                  mode='markers',
                                                  opacity=_OPACITY,
                                                  name='Outlier'))
            data.append(plotly.graph_objs.Scatter(x=dat['line1'][0],
                                                  y=dat['line1'][1],
                                                  text=dat['line1'][2],
                                                  mode='lines',
                                                  opacity=_OPACITY,
                                                  name='Helix'))
            data.append(plotly.graph_objs.Scatter(x=dat['line2'][0],
                                                  y=dat['line2'][1],
                                                  text=dat['line2'][2],
                                                  mode='lines',
                                                  opacity=_OPACITY,
                                                  name='Sheet'))

            layout = plotly.graph_objs.Layout(
                width=800,
                height=400,
                xaxis=dict(title='\u0394\u03B4CA - \u0394\u03B4CB'),
                yaxis=dict(title='\u0394\u03B4{}'.format(dat['atom'])),
                showlegend=True,
                hovermode='closest',
                title='{} [Offset = {} ppm]'.format(dat['atom'], dat['offset']))
            fig = plotly.graph_objs.Figure(data=data, layout=layout)
            plotly.offline.plot(fig, filename='LACS{}.html'.format(dat['atom']), auto_open=_AUTOOPEN)

if __name__ == "__main__":
    p = LACS_Viz()
    # p.read_lacs_output('data/2m07_1_A.lacs_out')

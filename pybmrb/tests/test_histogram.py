#!/usr/bin/env python3


from pybmrb import Histogram

histogram=Histogram()
def test_hist():
    data = histogram.hist(residue='CYS', show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 2
    assert len(data[0]) == len(data[1])


def test_hist2d():
    data = histogram.hist2d(residue='CYS', atom1='CB', atom2='CA', show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 2
    assert len(data[0]) == len(data[1])


def test_conditional_hist():
    data = histogram.conditional_hist(residue='THR', atom='N',
                                      filtering_rules=[('CB', 69.51), ('CA', 60.79), ('H', 8.13)],
                                      show_visualization=False)
    assert type(data) is tuple
    assert len(data) == 2
    assert len(data[0]) == len(data[1])
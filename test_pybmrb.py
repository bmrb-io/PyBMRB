#!/usr/bin/env python
from unittest import TestCase
import pybmrb as BV


# class TestHistogram(TestCase):
#     def test_get_histogram_api(self):
#         self.fail()
#
#     def test_get_conditional_histogram_api(self):
#         self.fail()
#
#     def test_get_histogram2d_api(self):
#         self.fail()
#
#     def test_single_2dhistogram(self):
#         self.fail()
#
#     def test_single_atom(self):
#         self.fail()
#
#     def test_multiple_atom(self):
#         self.fail()
#
#     def test_conditional_histogram(self):
#         self.fail()

class TestSpectra(TestCase):

    def setUp(self):
        self.Spectra = BV.Spectra()

    def test_get_entry(self):
        self.assertEqual(self.Spectra.get_entry(), [])
        self.assertGreater(len(self.Spectra.get_entry('15060')), 0)
        self.assertGreater(len(self.Spectra.get_entry(seq='MTFKLIINGKTLKGETTTEAVDAATAEKVFKQYFNDNGIDGEWTYDDATKTFTITE')),
                           0)
        self.assertGreater(len(self.Spectra.get_entry(bmrbid='15060',
                                                      seq='MTFKLIINGKTLKGETTTEAVDAATAEKVFKQYFNDNGIDGEWTYDDATKTFTITE')),
                           0)
        self.assertGreater(len(self.Spectra.get_entry(bmrbid='15060',
                                                      seq='MTFKLIINGKTLKGETTTEAVDAATAEKVFKQYFNDNGIDGEWTYDDATKTFTITE',
                                                      nn=10)),
                           0)

    # def test__load_pp_dict(self):
    #     self.fail()
    #
    # def test_predict_from_seq(self):
    #     self.fail()
    #
    # def test_n15_hsqc(self):
    #     self.fail()
    #
    # def test_plotn15_hsqc(self):
    #     self.fail()

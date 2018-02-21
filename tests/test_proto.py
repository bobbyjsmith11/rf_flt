#!/usr/bin/env python
"""
"""
from rf_flt import proto
import os
import random
from unittest import TestCase

class TestPrototypes(TestCase):
    """ Tests proto.py
    """

    def setUp(self):
        self.test_dir = os.path.dirname(os.path.abspath(__file__)) + "/"

    def test_butterworth_coefficients(self):
        # known good butterworth coefficients for 10 pole filter
        g_good_10 = [1.0,
                     0.3129,
                     0.9080,
                     1.4142,
                     1.7820,
                     1.9754,
                     1.9754,
                     1.7820,
                     1.4142,
                     0.9080,
                     0.3129,
                     1.0]
        g_test = proto.butterworth_coef(10)
        
        # choose a random coefficient
        n = int(random.random()*10) 
        self.assertAlmostEqual(g_good_10[n], round(g_test[n],4), places=3)
    
    def test_chebyshev_coefficients_0p5db_ripple(self):
        # known good chebyshev coefficients for 10 pole, 0.5dB ripple
        g_good_10 = [1.0,
                     1.7543,
                     1.2721,
                     2.6754,
                     1.3725,
                     2.7392,
                     1.3806,
                     2.7231,
                     1.3485,
                     2.5239,
                     0.8842,
                     1.9841]
        g_test = proto.chebyshev_coef(10, 0.5)
        
        # choose a random coefficient
        n = int(random.random()*10) 
        self.assertAlmostEqual(g_good_10[n], round(g_test[n],4), places=3)

    # def test_denormalize_proto(self):
    #     """
    #     """
    #     cmp_lst = proto.butterworth_coef(5)
    #     d = proto.denormalize_prototype(cmp_lst, 2e9, first_element='shunt')
    #     d_test = []
    #     for item in d:
    #         d_test.append(item['value'])
    #     d_good = [50.0,
    #               9.836316e-13,
    #               6.43795e-9,
    #               3.18309e-12,
    #               6.43795e-9,
    #               9.836316e-13,
    #               50.0]
    #     n = int(random.random()*7)
    #     self.assertAlmostEqual(d_test[n], d_good[n], places=3)
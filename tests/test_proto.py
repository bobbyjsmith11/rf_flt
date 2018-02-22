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

    # def test_connect_series_single_commponent(self):
    #     test_cmp = {'ref_des': 'X1', 'value': 1.0} 
    #     node = 1
    #     proto.connect_series(test_cmp, node, flt_type='lowpass')

    def test_connect_in_series_no_last_node(self):
        c1 = {'ref_des': 'X1', 'value': 1}
        c2 = {'ref_des': 'X2', 'value': 1}
        c3 = {'ref_des': 'X3', 'value': 1}
        cmps = [c1, c2, c3]
        d_good = {'X1': {'pins': [5, 6], 'value':1},
                  'X2': {'pins': [6, 7], 'value': 1},
                  'X3': {'pins': [7, 8], 'value': 1},
                  }
        d_test = proto.connect_in_series(cmps, 5)
        self.assertEqual(d_good, d_test)

    def test_connect_in_series_last_node_not_gnd(self):
        c1 = {'ref_des': 'X1', 'value': 1}
        c2 = {'ref_des': 'X2', 'value': 1}
        c3 = {'ref_des': 'X3', 'value': 1}
        cmps = [c1, c2, c3]
        d_good = {'X1': {'pins': [5, 9], 'value':1},
                  'X2': {'pins': [9, 7], 'value': 1},
                  'X3': {'pins': [7, 6], 'value': 1},
                  }
        d_test = proto.connect_in_series(cmps, 5, last_node=6)
        self.assertEqual(d_good, d_test)
    
    def test_connect_in_series_last_node_gnd(self):
        c1 = {'ref_des': 'X1', 'value': 1}
        c2 = {'ref_des': 'X2', 'value': 1}
        c3 = {'ref_des': 'X3', 'value': 1}
        cmps = [c1, c2, c3]
        d_good = {'X1': {'pins': [5, 6], 'value':1},
                  'X2': {'pins': [6, 7], 'value': 1},
                  'X3': {'pins': [7, 0], 'value': 1},
                  }
        d_test = proto.connect_in_series(cmps, 5, last_node=0)
        self.assertEqual(d_good, d_test)
    
    def test_connect_in_series_single_cmp(self):
        c1 = {'ref_des': 'X1', 'value': 1}
        d_good = {'X1': {'pins': [5, 6], 'value':1},
                  }
        d_test = proto.connect_in_series(c1, 5)
        self.assertEqual(d_good, d_test)
    
    def test_connect_in_parallel_no_last_node(self):
        c1 = {'ref_des': 'X1', 'value': 1}
        c2 = {'ref_des': 'X2', 'value': 1}
        c3 = {'ref_des': 'X3', 'value': 1}
        cmps = [c1, c2, c3]
        d_good = {'X1': {'pins': [5, 6], 'value':1},
                  'X2': {'pins': [5, 6], 'value': 1},
                  'X3': {'pins': [5, 6], 'value': 1},
                  }
        d_test = proto.connect_in_parallel(cmps, 5)
        self.assertEqual(d_good, d_test)
    
    def test_connect_in_parallel_last_node_not_gnd(self):
        c1 = {'ref_des': 'X1', 'value': 1}
        c2 = {'ref_des': 'X2', 'value': 1}
        c3 = {'ref_des': 'X3', 'value': 1}
        cmps = [c1, c2, c3]
        d_good = {'X1': {'pins': [5, 10], 'value':1},
                  'X2': {'pins': [5, 10], 'value': 1},
                  'X3': {'pins': [5, 10], 'value': 1},
                  }
        d_test = proto.connect_in_parallel(cmps, 5, last_node=10)
        self.assertEqual(d_good, d_test)

    def test_connect_in_parallel_last_node_gnd(self):
        c1 = {'ref_des': 'X1', 'value': 1}
        c2 = {'ref_des': 'X2', 'value': 1}
        c3 = {'ref_des': 'X3', 'value': 1}
        cmps = [c1, c2, c3]
        d_good = {'X1': {'pins': [5, 0], 'value':1},
                  'X2': {'pins': [5, 0], 'value': 1},
                  'X3': {'pins': [5, 0], 'value': 1},
                  }
        d_test = proto.connect_in_parallel(cmps, 5, last_node=0)
        self.assertEqual(d_good, d_test)
    





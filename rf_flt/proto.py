#!/usr/bin/env python
"""
===========================
proto.py
===========================
:Description:
    Module for calculating the various prototype component values, i.e.
    Butterworth, Chebyshev, etc.

"""
import math

def butterworth_coef(poles, terms=2):
    """
    Return a list of Butterworth (a.k.a Maximally Flat) coefficients,
    given number of poles.
    These are the normalized component values (Z0=1 Ohm and omega_c=1 radian/second)
    """
    if terms == 2:
        g=[]
        g.append(1.0)
        for k in range(1,poles+1):
            g.append(2*math.sin((2*k-1)*math.pi/(2*poles)))
        g.append(1.0)
    else:
        g = "this algorithm needs work"
    return g

def chebyshev_coef(poles, ripple_db):
    """
    Return a list of Chebyshev (a.k.a. Equal Ripple) coefficients, 
    given number of poles and  passband ripple. 
    These are the normalized component values (Z0=1 Ohm and omega_c=1 radian/second)
    """ 
    beta = math.log(1/math.tanh(ripple_db/17.37))  
    gamma = math.sinh(beta/(2*poles))   
    A=[]
    B=[]
    A.append(0)
    B.append(0)
    g=[]
    for k in range(1,poles+1):
        A.append(math.sin((((2*k)-1)*math.pi)/(2*poles)))
        B.append(gamma**2 + (math.sin(k*math.pi/poles))**2)

    g.append(1.0)
    g.append(2*A[1]/gamma)

    for k in range(2,poles+1):
        g.append((4.0*A[k]*A[k-1])/(B[k-1]*g[k-1]))
    if poles%2==0:
        g.append((1/(math.tanh(beta/4)))**2)
    else:
        g.append(1.0)
    return g

def denormalize_prototype(cmp_lst, fc, R0=50, first_element='shunt'):
    """ return the denormalized list of components
    :Parameters:
        - cmp_lst (list) - normalized list of components (coefficients)
        - fc (int or float) - cutoff frequency in Hz
        - R0 (impedance) - network impedance, default to 50 Ohms
    :Returns:
        list of dicts
    """
    ret = {}
    
    # source terminatin resistance
    next_node = 1
    # ret.append({'ref_des': 'RS',
    #             'value': cmp_lst[0]*R0,
    #             'pins': ["net_" + str(next_node), 'GND']})
    ret['RS'] = {'value': cmp_lst[0]*R0,
                'pins': ["net_" + str(next_node), 'GND']}
    for i in range(1, len(cmp_lst)-1):
        if first_element == "shunt":
            if (i % 2):
                ref_des = "C" + str(i)
                val = cmp_lst[i]/(2*math.pi*fc*R0)
                plus_node = 'net_' + str(next_node)
                minus_node = 'GND '
                # it's odd, shunt C
            else:
                # it's even, series L
                ref_des = "L" + str(i)
                val = cmp_lst[i]*R0/(2*math.pi*fc)
                plus_node = 'net_' + str(next_node)
                next_node += 1
                minus_node = 'net_' + str(next_node)
        else:
            if (i % 2):
                ref_des = "L" + str(i)
                val = cmp_lst[i]*R0/(2*math.pi*fc)
                plus_node = 'net_' + str(next_node)
                next_node += 1
                minus_node = 'net_' + str(next_node)
            else:
                ref_des = "C" + str(i)
                val = cmp_lst[i]/(2*math.pi*fc*R0)
                plus_node = 'net_' + str(next_node)
                minus_node = 'GND' 
        ret[ref_des] = {'value': val,
                        'pins': [str(plus_node), str(minus_node)]}
    
    ret['RL'] = {'value': cmp_lst[-1]*R0,
                'pins': ["net_" + str(next_node), 'GND']}
    return ret



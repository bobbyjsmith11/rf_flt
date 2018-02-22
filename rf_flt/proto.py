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

class Filter(object):
    """
    """
    def __init__(self, num_poles, type='butter'):
        """
        Parameters
            num_poles (int) - number of elements in the filter
            type (str) - 'butter', 'cheby'

        """
        pass

        
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

def lowpass_transform(g, fc, R0=50, first_element='shunt'):
    """ return the denormalized list of components
    :Parameters:
        - g (list) - normalized list of components (coefficients)
        - fc (int or float) - cutoff frequency in Hz
        - R0 (impedance) - network impedance, default to 50 Ohms
    :Returns:
        list of dicts
    """
    ret = {}
    
    # source terminating resistance
    next_node = 1
    # ret.append({'ref_des': 'RS',
    #             'value': g[0]*R0,
    #             'pins': ["net_" + str(next_node), 'GND']})
    ret['RS'] = {'value': g[0]*R0,
                'pins': ["net_" + str(next_node), 'GND']}
    for i in range(1, len(g)-1):
        if first_element == "shunt":
            if (i % 2):
                ref_des = "C" + str(i)
                val = g[i]/(2*math.pi*fc*R0)
                plus_node = 'net_' + str(next_node)
                minus_node = 'GND '
                # it's odd, shunt C
            else:
                # it's even, series L
                ref_des = "L" + str(i)
                val = g[i]*R0/(2*math.pi*fc)
                plus_node = 'net_' + str(next_node)
                next_node += 1
                minus_node = 'net_' + str(next_node)
        else:
            if (i % 2):
                ref_des = "L" + str(i)
                val = g[i]*R0/(2*math.pi*fc)
                plus_node = 'net_' + str(next_node)
                next_node += 1
                minus_node = 'net_' + str(next_node)
            else:
                ref_des = "C" + str(i)
                val = g[i]/(2*math.pi*fc*R0)
                plus_node = 'net_' + str(next_node)
                minus_node = 'GND' 
        ret[ref_des] = {'value': val,
                        'pins': [str(plus_node), str(minus_node)]}
    
    ret['RL'] = {'value': g[-1]*R0,
                'pins': ["net_" + str(next_node), 'GND']}
    return ret

def transform_prototype(gs, 
                        fc, 
                        R0=50, 
                        bw_ratio=0.1,
                        first_element='shunt', 
                        flt_type='lowpass'):
    """
    Parameters
        gs (list of floats) - prototype filter coefficients
        fc (float) - cutoff or center frequency (depending on flt_type)
        R0 (float) - filter impedance (default = 50)
        bw_ratio (float) - relative bandwidth (0.1 = 10%).
                            This only applies to bandpass or bandreject
        first_element (str) - <'series', 'shunt'>
        flt_type (str) - <'lowpass', 'highpass', 'bandpass', 'bandreject'>
    """
    ret = {}
    
    cmp_lst = [] 
    cmp_lst.append("RS")
    for i in range(1, len(gs)-1):
        if first_element == "shunt":
            if (i % 2):
                # it's odd, shunt C
                cmp = transform_shunt_coeff(gs[i], 
                                           fc, 
                                           R0=R0, 
                                           bw_ratio=bw_ratio, 
                                           flt_type=flt_type)
            else:
                cmp = transform_series_coeff(gs[i], 
                                            fc, 
                                            R0=R0, 
                                            bw_ratio=bw_ratio, 
                                            flt_type=flt_type)
        else:
            if (i % 2):
                cmp = transform_series_coeff(gs[i], 
                                            fc, 
                                            R0=R0, 
                                            bw_ratio=bw_ratio, 
                                            flt_type=flt_type)
            else:
                cmp = transform_shunt_coeff(gs[i], 
                                           fc, 
                                           R0=R0, 
                                           bw_ratio=bw_ratio, 
                                           flt_type=flt_type)
        for item in cmp:
            item['ref_des'] = item['ref_des'] + str(i) 
        cmp_lst.append(cmp)
        
    cmp_lst.append("RL")
    return cmp_lst

def connect_series(cmp_list, next_node, flt_type='lowpass'):
    """
    returns tuple(list of dicts, next_node) 
    """
    d_ret = {}
    if flt_type == 'lowpass' or flt_type == 'highpass':
        ref_des = cmp_list[0]['ref_des']
        val = cmp_list[0]['value']
        plus_node = 'net_' + str(next_node)
        next_node += 1
        minus_node = 'net_' + str(next_node)
        d_ret[ref_des] = {'value': val,
                        'pins': [str(plus_node), str(minus_node)]}
    elif flt_type == 'bandpass':
        for cmp in cmp_list:
            ref_des = cmp['ref_des']
            val = cmp['value']
            plus_node = 'net_' + str(next_node)
            next_node += 1
            minus_node = 'net_' + str(next_node)
            d_ret[ref_des] = {'value': val,
                            'pins': [str(plus_node), str(minus_node)]}
    elif flt_type == 'bandreject':
        for cmp in cmp_list:
            ref_des = cmp['ref_des']
            val = cmp['value']
            plus_node = 'net_' + str(next_node)
            minus_node = 'net_' + str(next_node + 1)
            d_ret[ref_des] = {'value': val,
                            'pins': [str(plus_node), str(minus_node)]}
        next_node += 1

    return (d_ret, next_node)

def transform_series_coeff(g, fc, R0=50, bw_ratio=0.1, flt_type='lowpass'):
    """ 
    perform the appropriate series transform on the given normalized 
    filter coefficient. Return  a list of tuples ('L' or 'C' as str, value of component as float)
    """
    lst = []
    if flt_type == 'lowpass':
        val = R0*g/(2*math.pi*fc)
        lst.append({'ref_des':'L', 'value':val}) 

    elif flt_type == 'highpass':                
        val = 1.0/(g*R0*2*math.pi*fc) 
        lst.append({'ref_des':'C', 'value':val}) 

    elif flt_type == 'bandpass':
        val = R0*g/(2*math.pi*fc*bw_ratio)
        lst.append({'ref_des':'L', 'value':val}) 
        val = bw_ratio/(R0*2*math.pi*fc*g)
        lst.append({'ref_des':'C', 'value':val}) 
    
    elif flt_type == 'bandreject':
        val = R0*g*bw_ratio/(2*math.pi*fc)
        lst.append({'ref_des':'L', 'value':val}) 
        val = 1.0/(R0*g*2*math.pi*fc*bw_ratio)
        lst.append({'ref_des':'C', 'value':val}) 
    else:
        raise ValueError("not a valid filter type")
    return lst

def transform_shunt_coeff(g, fc, R0=50, bw_ratio=0.1, flt_type='lowpass'):
    """
    perform the appropriate shunt transform on the given normalized 
    filter coefficient. Return  a list of tuples ('L' or 'C' as str, value of component as float)
    """
    lst = []
    if flt_type == 'lowpass':
        val = g/(R0*2*math.pi*fc)
        lst.append({'ref_des':'C', 'value':val}) 
    
    elif flt_type == 'highpass':
        val = g*R0/(2*math.pi*fc)
        lst.append({'ref_des':'L', 'value':val}) 

    elif flt_type == 'bandpass':
        val = bw_ratio*R0/(2*math.pi*fc*g)
        lst.append({'ref_des':'L', 'value':val}) 
        val = g/(2*math.pi*fc*bw_ratio*R0) 
        lst.append({'ref_des':'C', 'value':val}) 

    elif flt_type == 'bandreject':
        val = R0/(g*2*math.pi*fc*bw_ratio)
        lst.append({'ref_des':'L', 'value':val}) 
        val = g*bw_ratio/(R0*2*math.pi*fc)
        lst.append({'ref_des':'C', 'value':val}) 
    else:
        raise ValueError("not a valid filter type")
    return lst
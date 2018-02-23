#!/usr/bin/env python
"""
===========================
proto.py
===========================
:Description:
    Module for calculating the various prototype component values, i.e.
    Butterworth, Chebyshev, etc.
:Example Usage:

    >>> from rf_flt import proto
    >>> flt = proto.Filter(5, type='butter')
    >>> d = flt.get_bandpass(2e9, 0.1)
    >>> from rf_flt import spice_net
    >>> net = spice_net.SpiceNetwork(d)
    >>> net.write_netlist(fo="my_netlist.cir", title="MyBPF")

"""
import math
from . import spice_net

class Filter(object):
    """
    """
    def __init__(self, 
                 num_poles, 
                 type='butter', 
                 ripple=0.5,
                 R0=50):
        """
        Parameters
            num_poles (int) - number of poles in filter
            type (str) - <'butter', 'cheby'>
            ripple (float) - for Chebyshev only this specifies the ripple in dB
            R0 (float) - impedance in Ohms
        """
        self.type = type
        self.R0 = R0
        self.num_poles = num_poles
        self.ripple = ripple

    @property
    def gs(self):
        return self.get_coefficients()

    def get_coefficients(self):
        if self.type == 'butter': 
            g = butterworth_coef(self.num_poles)
        elif self.type == 'cheby':
            g = chebyshev_coef(self.num_poles, self.ripple)
        return g

    def get_lowpass(self, fc=1.0, first_element='shunt', normalized=False):
        """ 
        Return a spice_net.SpiceNetwork object for the given filter
        :Parameters:
            fc (float) - cutoff frequency in Hz
            first_element (str) - <'series'|'shunt'>
            normalized (bool) - set True to get values normalized to R0=1 Ohm and fc=1 rad/sec
        """ 
        if not normalized:
            cmps = transform_prototype(self.gs, 
                                       fc, 
                                       R0=self.R0, 
                                       first_element=first_element, 
                                       flt_type='lowpass')
        else:
            cmps = transform_prototype_normalized(self.gs, 
                                                  first_element=first_element, 
                                                  flt_type='lowpass')
        cmp_dict = connect_highpass_lowpass(cmps, first_element=first_element)
        net = spice_net.SpiceNetwork(cmp_dict)
        return net

    def get_highpass(self, fc=1.0, first_element='shunt', normalized=False):
        """ 
        Return a spice_net.SpiceNetwork object for the given filter
        :Parameters:
            fc (float) - cutoff frequency in Hz
            first_element (str) - <'series'|'shunt'>
            normalized (bool) - set True to get values normalized to R0=1 Ohm and fc=1 rad/sec
        """ 
        if not normalized:
            cmps = transform_prototype(self.gs, 
                                       fc, 
                                       R0=self.R0, 
                                       first_element=first_element, 
                                       flt_type='highpass')
        else:
            cmps = transform_prototype_normalized(self.gs, 
                                                  first_element=first_element, 
                                                  flt_type='highpass')
        cmp_dict = connect_highpass_lowpass(cmps, first_element=first_element)
        net = spice_net.SpiceNetwork(cmp_dict)
        return net

    def get_bandpass(self, fc=1.0, bw_ratio=0.2, first_element='shunt', normalized=False):
        """ 
        Return a spice_net.SpiceNetwork object for the given filter
        :Parameters:
            fc (float) - center frequency in Hz
            bw_ratio (float) - bandwidth (0.1=10%)
            first_element (str) - <'series'|'shunt'>
            normalized (bool) - set True to get values normalized to R0=1 Ohm and fc=1 rad/sec
        """ 
        if not normalized:
            cmps = transform_prototype(self.gs, 
                                       fc,
                                       bw_ratio=bw_ratio, 
                                       R0=self.R0, 
                                       first_element=first_element, 
                                       flt_type='bandpass')
        else:
            cmps = transform_prototype_normalized(self.gs,
                                                  bw_ratio=bw_ratio, 
                                                  first_element=first_element, 
                                                  flt_type='bandpass')
        cmp_dict = connect_bandpass(cmps, first_element=first_element)
        net = spice_net.SpiceNetwork(cmp_dict)
        return net
        

    def get_bandstop(self, fc=1.0, bw_ratio=0.2, first_element='shunt', normalized=False):
        """ 
        Return a spice_net.SpiceNetwork object for the given filter
        :Parameters:
            fc (float) - center frequency in Hz
            bw_ratio (float) - bandwidth (0.1=10%)
            first_element (str) - <'series'|'shunt'>
            normalized (bool) - set True to get values normalized to R0=1 Ohm and fc=1 rad/sec
        """ 
        if not normalized:
            cmps = transform_prototype(self.gs, 
                                       fc,
                                       bw_ratio=bw_ratio, 
                                       R0=self.R0, 
                                       first_element=first_element, 
                                       flt_type='bandstop')
            cmps = transform_prototype_normalized(self.gs,
                                                  bw_ratio=bw_ratio, 
                                                  first_element=first_element, 
                                                  flt_type='bandstop')
        cmp_dict = connect_bandstop(cmps, first_element=first_element)
        net = spice_net.SpiceNetwork(cmp_dict)
        return net

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

def transform_prototype_normalized(gs, bw_ratio=0.1, first_element='shunt', flt_type='lowpass'):
    """
    """
    fc = 1.0/(2*math.pi)
    cmp_list = transform_prototype(gs, fc, R0=1, bw_ratio=bw_ratio, first_element=first_element, flt_type=flt_type)
    return cmp_list

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
                            This only applies to bandpass or bandstop
        first_element (str) - <'series', 'shunt'>
        flt_type (str) - <'lowpass', 'highpass', 'bandpass', 'bandstop'>
    """
    cmp_lst = [] 
    for i in range(1, len(gs)-1):
        if first_element == "shunt":
            if (i % 2):
                # it's odd, shunt 
                cmp = transform_shunt_coeff(gs[i], 
                                           fc, 
                                           R0=R0, 
                                           bw_ratio=bw_ratio, 
                                           flt_type=flt_type)
            else:
                # it's even, series 
                cmp = transform_series_coeff(gs[i], 
                                            fc, 
                                            R0=R0, 
                                            bw_ratio=bw_ratio, 
                                            flt_type=flt_type)
        else:
            if (i % 2):
                # it's odd, series 
                cmp = transform_series_coeff(gs[i], 
                                            fc, 
                                            R0=R0, 
                                            bw_ratio=bw_ratio, 
                                            flt_type=flt_type)
            else:
                # it's even, shunt 
                cmp = transform_shunt_coeff(gs[i], 
                                           fc, 
                                           R0=R0, 
                                           bw_ratio=bw_ratio, 
                                           flt_type=flt_type)
        for item in cmp:
            item['ref_des'] = item['ref_des'] + str(i) 
        cmp_lst.append(cmp)
        
    return cmp_lst

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
        val = 1.0/(R0*2*math.pi*fc*g) 
        lst.append({'ref_des':'C', 'value':val}) 

    elif flt_type == 'bandpass':
        val = R0*g/(2*math.pi*fc*bw_ratio)
        lst.append({'ref_des':'L', 'value':val}) 
        val = bw_ratio/(R0*2*math.pi*fc*g)
        lst.append({'ref_des':'C', 'value':val}) 
    
    elif flt_type == 'bandstop':
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
        val = R0/(2*math.pi*fc*g)
        lst.append({'ref_des':'L', 'value':val}) 

    elif flt_type == 'bandpass':
        val = bw_ratio*R0/(2*math.pi*fc*g)
        lst.append({'ref_des':'L', 'value':val}) 
        val = g/(2*math.pi*fc*bw_ratio*R0) 
        lst.append({'ref_des':'C', 'value':val}) 

    elif flt_type == 'bandstop':
        val = R0/(g*2*math.pi*fc*bw_ratio)
        lst.append({'ref_des':'L', 'value':val}) 
        val = g*bw_ratio/(R0*2*math.pi*fc)
        lst.append({'ref_des':'C', 'value':val}) 
    else:
        raise ValueError("not a valid filter type")
    return lst

def connect_in_series(cmp_list, first_node, last_node=None):
    """
    Parameters
        cmp_list (list of dicts) - each dict has keys: 'ref_des', 'value'
    Returns dict with keys of 'ref_des' with values of dicts
    with the following keys: 
        'pins' : list of net names
        'value' : float
    """
    if not isinstance(cmp_list, list):
        cmp_list = [cmp_list]
    nodes = []
    for i in range(len(cmp_list) + 1):                  # create a list of node numbers but missing the last node
        nodes.append(first_node + i)
    if last_node != None:                               # check if user wants to name the last node
        if last_node in nodes:                          # if last_node is in nodes, replace it with next available node number
            idx = nodes.index(last_node)
            nodes[idx] = max(nodes) + 1
        nodes[-1] = last_node

    d_ret = {}
    for i in range(len(nodes) - 1):
        ref_des = cmp_list[i]['ref_des']
        val = cmp_list[i]['value']
        p_node = nodes[i]
        m_node = nodes[i+1]
        d_ret[ref_des] = {'value': val,
                        'pins': [p_node, m_node]}
    return d_ret

def connect_in_parallel(cmp_list, first_node, last_node=None):
    """
    Parameters
        cmp_list (list of dicts) - each dict has keys: 'ref_des', 'value'
    Returns dict with keys of 'ref_des' with values of dicts
    with the following keys: 
        'pins' : list of net names
        'value' : float
    """
    if not isinstance(cmp_list, list):
        cmp_list = [cmp_list]
    d_ret = {}
    for cmp in cmp_list:
        ref_des = cmp['ref_des']
        val = cmp['value']
        if last_node == None:
           last_node = first_node + 1 
        d_ret[ref_des] = {'value': val,
                        'pins': [first_node, last_node]}
    return d_ret

def connect_highpass_lowpass(cmp_list, first_element='shunt'):
    """ return a dict:
            keys are ref_des
            values are {'pins': [net1, net2], 'value': part value}
        The nodes are connected for a lowpass or highpass
    """
    d_ret = {}
    current_node = 1
    for i in range(len(cmp_list)):
        if first_element == 'shunt':
            if (i % 2):
                # it's series in series
                d = connect_in_series(cmp_list[i], current_node)
                current_node = current_node + 1
            else:
                # it's series to ground
                d = connect_in_series(cmp_list[i], current_node, last_node=0)
        else:
            if (i % 2):
                # it's series to ground
                d = connect_in_series(cmp_list[i], current_node, last_node=0)
            else:
                # it's series in series
                d = connect_in_series(cmp_list[i], current_node)
                current_node = current_node + 1
        d_ret.update(d)
    return d_ret

def connect_bandpass(cmp_list, first_element='shunt'):
    """ return a dict:
            keys are ref_des
            values are {'pins': [net1, net2], 'value': part value}
        The nodes are connected for a bandpass configuration
    """
    d_ret = {}
    current_node = 1
    for i in range(len(cmp_list)):
        if first_element == 'shunt':
            if (i % 2):
                # it's series in series
                d = connect_in_series(cmp_list[i], current_node)
                current_node = current_node + 2 
            else:
                # it's parallel to ground
                d = connect_in_parallel(cmp_list[i], current_node, last_node=0)
        else:
            if (i % 2):
                # it's parallel to ground
                d = connect_in_parallel(cmp_list[i], current_node, last_node=0)
            else:
                # it's series in series
                d = connect_in_series(cmp_list[i], current_node)
                current_node = current_node + 2
        d_ret.update(d)
    return d_ret

def connect_bandstop(cmp_list, first_element='shunt'):
    """ return a dict:
            keys are ref_des
            values are {'pins': [net1, net2], 'value': part value}
        The nodes are connected for a bandpass configuration
    """
    d_ret = {}
    current_node = 1
    for i in range(len(cmp_list)):
        if first_element == 'shunt':
            if (i % 2):
                # it's parallel in series
                d = connect_in_parallel(cmp_list[i], current_node, last_node=current_node+2)
                current_node = current_node + 2
            else:
                # it's series to ground
                d = connect_in_series(cmp_list[i], current_node, last_node=0)
        else:
            if (i % 2):
                # it's series to ground
                d = connect_in_series(cmp_list[i], current_node, last_node=0)
            else:
                # it's parallel in series
                d = connect_in_parallel(cmp_list[i], current_node, last_node=current_node+2)
                current_node = current_node + 2
        d_ret.update(d)
    return d_ret


def richards_transform(cmp_list, fc=1.0, normalized=True):
    """ 
    performa Richards Transormation on a component dictionary and
    return new component dictionary with only transmission lines 
    """
    d_ret = {}
    pin_list = []
    for ref_des in cmp_list:
        idx = ref_des[1:]
        val = cmp_list[ref_des]['value']
        p1 = cmp_list[ref_des]['pins'][0]
        pin_list.append(p1)
        p2 = cmp_list[ref_des]['pins'][1]
        pin_list.append(p2)
        if ref_des.startswith("L"):
            val_str = 'Z0=' + str(val) + ' F=' + str(fc) + ' NL=0.125'
            d_ret['T'+str(idx)] = {'pins': [p1, 'SC', p2, 'SC'],
                               'value': val_str}            
        elif ref_des.startswith("C"):
            val_str = 'Z0=' + str(1/val) + ' F=' + str(fc) + ' NL=0.125'
            d_ret['T'+str(idx)] = {'pins': [p1, 'OC', p2, 'OC'],
                               'value': val_str}            
        else:
            raise ValueError(ref_des + " not supported in Richards Transform")
    next_pin = max(pin_list) + 1
    for ref in d_ret:
        if d_ret[ref]['pins'][1] == 'SC':
            d_ret[ref]['pins'][1] = next_pin 
            d_ret[ref]['pins'][3] = next_pin 
            next_pin += 1
        elif d_ret[ref]['pins'][1] == 'OC':
            d_ret[ref]['pins'][1] = next_pin 
            next_pin += 1
            d_ret[ref]['pins'][3] = next_pin 
            next_pin += 1
    return d_ret


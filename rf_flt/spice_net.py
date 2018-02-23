#!/usr/bin/env python
"""
=================================
spice_net.py
=================================
:Description:
    Captures filter networks in netlist format similar to SPICE

Two data structures are widely used in this module

    1. cmp_list 
       This is a list of dictionaries. Each dictionary has the following key-value paris
            'ref_des': (str) reference designator (R5 for example)
            'value': (float) value of the element (1.5e-9 could mean 1.5nH or 1.5nF depending on the type of element) 
            'pins': (list of str) - the net name that each pin is connected to. 
            >>> cmp_list[1]['ref_des'] # returns
            >>> 'C1'
            This means the the second element in the list is for C1, and if we type:
            >>> cmp_list[1]['pins'][1] # returns
            >>> 'GND'
            This means that the second element in the list (C1) is connected to GND

    2. d_map
        This is a dictionary. The keys are the net names and each value is a list of 
        components with the pin number that is attached to this net. For example,
        >>> d_map['net_1']
        >>> ['RS.1', 'C1.1', 'L2.1']
        This means that net_1 is connected to RS pin 1, C1 pin 1 and L2 pin 1

"""

class SpiceNetwork(object):
    """
    """
    def __init__(self, parts_dict):
        self.components = parts_dict

    def write_netlist(self, fo=None, title="MYCIRCUIT.CIR"):
        if fo:
            if isinstance(fo, str):
                fo = open(fo, 'w')

        title_line = title
        print(title_line)
        if fo:
            fo.write(title_line + "\n")
        for ref_des in self.components:
            dat = []
            dat.append(ref_des)
            for node in self.components[ref_des]['pins']:
                dat.append(str(node))
            dat.append(self.components[ref_des]['value'])
            line = ("".join("{:<8}".format(x) for x in dat))
            print(line)
            if fo:
                fo.write(line + "\n")
        print(".END")
        if fo:
            fo.write('.END\n')

def read_netlist(fi):
    """ 
    Parameters
        fi (str) - input netlist file
    """
    with open(fi, 'r') as f:
        lines = f.readlines()
    lines = [x.strip() for x in lines]
    cmp_dict = {}
    for line in lines:
        if line.startswith("*"):    
            # it's a comment
            continue
        elif line.startswith(".end") or line.startswith(".END"):
            break
        else:
            ref_des, value, pins = parse_comp_line(line)
            cmp_dict[ref_des] = {'value': value, 'pins': pins}
    return cmp_dict

def parse_comp_line(cmp_str):
    """
    """
    cmp_ar = cmp_str.split()
    ref_des = cmp_ar.pop(0)
    value = cmp_ar.pop()
    pins = []
    for pin in cmp_ar:
        pins.append(pin)
    return ref_des, value, pins 

def get_netlist_comp(cmp_list):
    """
    Keyword Arguments
        cmp_list (list) - list of components as dicts
    
    Returns
         dict wit the net names as keys. the values of each key is a list of component ref_des with pin number
    """
    nets = []
    d = {}
    for cmp in cmp_list:
        for i in range(len(cmp['pins'])):
            if cmp['pins'][i] not in d.keys():
                d[cmp['pins'][i]] = []
            d[cmp['pins'][i]].append(cmp['ref_des'] + "." + str(i+1))
            nets.append(cmp['pins'][i])
    nets = list(set(nets))

    return d

def find_components_with_net(cmp_list, net="GND"):
    """
    find all components connected to the given net name. 
    return list
    Parameters
        cmp_list (list) - list of components as dicts
    """
    d_map = get_netlist_comp(cmp_list)
    ref_des_lst = []
    for cmp_pin in d_map[net]:
        ref_des_lst.append(cmp_pin.split(".")[0])
    return ref_des_lst

def find_parallel_components(cmp_list, ref_des):
    """
    find all components whose pins are connected to the same nets
    Parameters
        cmp_list (list) - list of components as dicts
    """
    for cmp in cmp_list:
        pass
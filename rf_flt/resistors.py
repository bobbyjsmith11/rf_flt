#!/usr/bin/env python
"""
===========================
resistors.py
===========================
:Description:
    various calculations involving resistor networks

"""
import numpy as np

def r_parallel(r1, r2):
    """
    return the parallel combo of two resistors
    ____________________
        |       |
        /       /
        \       \
        /       /
     r1 \    r2 \    rp = r1 || r2
        /       /
        \       \
        |       |
       ---     ---
        -       -
    """
    return r1*r2/(r1+r2)

def find_r1_parallel(r2, rp):
    """
    return the value of r1 when combined with
    r2 gives the parallel resistance of rp
    ____________________
        |       |
        /       /
        \       \
        /       /
     r1 \    r2 \    rp = r1 || r2
        /       /
        \       \
        |       |
       ---     ---
        -       -
    """
    if rp > r2:
        raise ValueError("rp cannot be greater than r2")
    return rp/(1 - rp/r2)


def get_r1r2_ratio(vcc, vout):
    """
    return the resistor ratio r1/r2
    given the supply voltage vcc and the
    desired vout
    _________ vcc
        |    
        /    
        \    
        /    
     r1 \    
        /    
        \    
        |    
        | --- vout 
        /
        \
        /
     r2 \ 
        /
        \
        |
       ---
        -
    """
    return (vcc-vout)/vout


def get_r2_thevenin(r1, rp, r1r2):
    """
    """
    return (rp + r1r2*rp)/r1r2

def get_pi_atten(a_db, z0=50):
    """
    return r1, r2 values for given attenuation
    :args:
        :a_db (float): attenuation in db
        :z0 (float): impedance
    :returns:
        (r1, r2) where r1 is shunt and r2 is series
    """
    r1 = z0*((10**(a_db/20)+1)/(10**(a_db/20)-1))
    r2 = (z0/2)*(10**(a_db/20) - 1/(10**(a_db/20)))
    return r1, r2


def get_tee_atten(a_db, z0=50):
    """
    return r1, r2 values for given attenuation
    :args:
        :a_db (float): attenuation in db
        :z0 (float): impedance
    """
    r1 = z0*((10**(a_db/20)-1)/(10**(a_db/20)+1))
    r2 = 2*z0*10**(a_db/20)/(10**(a_db/10) - 1)
    return r1, r2


def get_z_pi(r1, r2, z0=50):
    """
    return the impedance of the attenuator given r1 and r2
    """
    return r_parallel(r1, r2 + r_parallel(r1, z0))

def get_z_tee(r1, r2, z0=50):
    """
    return the impedance of the attenuator given r1 and r2
    """
    return r1 + r_parallel(r2, r1+ z0)


def calc_pi_attenuation(r1, r2, z0=50):
    """
    :Args:
        :r1 (float): shunt resistor
        :r2 (float): series resistor
    """
    r1p = (r1* z0) / (r1 + z0)
    zL = (r1 * (r1p + r2)) / (r1 + (r1p + r2))
    rL = ((zL - z0) / (zL + z0))
    rLdB = 20 * np.log10(np.abs(rL))
    Pin = 1
    vsource = 2 * np.sqrt(z0 * Pin)
    Iin = vsource / (z0 + zL)
    Vin = vsource - (Iin * z0)
    Vou = Vin * (r1p / (r1p + r2))
    Pou = Vou * Vou / z0
    gain = 10 * np.log10(Pou / Pin)
    atten = -1 * gain
    return atten, rLdB

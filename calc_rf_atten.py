"""

"""

import argparse
from rf_flt import resistors


def get_atten(atten, atten_type):
    pass

def show_pi(atten, z0):
    """
    """
    print( """
              R2 
    _________|_|_|________
        |             |
        -             - 
     R1 -         R1  -
        -             -
        |             |
       ---           ---
        -             -
    """)
    r1, r2 = resistors.get_pi_atten(a_db=atten, z0=z0)
    print("R1: {:<.2f}".format(r1))
    print("R2: {:<.2f}".format(r2))

def show_tee(atten, z0):
    """
    """
    print( """
       R1              R1 
    __|_|_|__________|_|_|___
               |
               - 
           R2  -
               -
               |
              ---
               -
    """)
    r1, r2 = resistors.get_tee_atten(a_db=atten, z0=z0)
    print("R1: {:<.2f}".format(r1))
    print("R2: {:<.2f}".format(r2))
if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("atten", type=float, help="attenuation in dB")
    parser.add_argument("--type", type=str, help="pi or tee. pi is default", default='pi')
    parser.add_argument("--z0", type=float, help="impedance", default=50)
    args = parser.parse_args()
    if args.type == 'pi':
        show_pi(args.atten, args.z0)
    else:
        show_tee(args.atten, args.z0)

"""

"""
import argparse
import numpy as np
import pandas as pd
from tabulate import tabulate

class Block(object):
    def __init__(self, g_db=0, nf_db=0,
                 oip3_dbm=999, oip2_dbm=999,
                 op1db_dbm=999, name='stage'):
        """

        """
        self.g_db = float(g_db)
        self.nf_db = float(nf_db)
        self.oip3_dbm = float(oip3_dbm)
        self.oip2_dbm = float(oip2_dbm)
        self.op1db_dbm = float(op1db_dbm)
        self.name = name


    def get_dict(self):
        """
        """
        return {
            'g_db': self.g_db,
            'nf_db': self.nf_db,
            'oip3_dbm': self.oip3_dbm,
            'iip3_dbm': self.iip3_dbm,
            'oip2_dbm': self.oip2_dbm,
            'iip2_dbm': self.iip2_dbm,
            'op1db_dbm': self.op1db_dbm,
            'ip1db_dbm': self.ip1db_dbm,
        }


    def __str__(self):
        return "Block(name={}, g_db={}, nf_db={}, oip3_dbm={}, oip2_dbm={}, op1db_dbm={})".format(self.name,
                                                                                                  self.g_db,
                                                                                                  self.nf_db,
                                                                                                  self.oip3_dbm,
                                                                                                  self.oip2_dbm,
                                                                                                  self.op1db_dbm)


    @property
    def iip2_dbm(self):
        return self.oip2_dbm - self.g_db
    @iip2_dbm.setter
    def iip2_dbm(self, val):
        self.oip2_dbm = val + self.g_db

    @property
    def iip3_dbm(self):
        return self.oip3_dbm - self.g_db
    @iip3_dbm.setter
    def iip3_dbm(self, val):
        self.oip3_dbm = val + self.g_db

    @property
    def ip1db_dbm(self):
        return self.op1db_dbm - self.g_db
    @ip1db_dbm.setter
    def ip1db_dbm(self, val):
        self.op1db_dbm = val + self.g_db

########################################
# CONVERSIONS
########################################
def db_to_ratio(db_in):
    """
    return the linear ratio given the db
    """
    return 10**(db_in/10)


def ratio_to_db(ratio):
    """
    return db of the given ratio
    """
    return 10*np.log10(ratio)


def dbm_to_w(dbm):
    """
    return watts from dBm
    """
    return (10**(dbm/10))*0.001


def w_to_dbm(w):
    """
    return dBm from watts
    """
    return 10*np.log10(w/0.001)
########################################
# END OF CONVERSIONS
########################################


def get_cascade(blks):
    """
    return a list of noise figures
    """
    d = {
        'gain_c': [],
        'nf_c': [],
        'oip3_c': [],
        'oip2_c': [],
        'op1db_c': []
    }
    g_running = blks[0].g_db
    d['gain_c'].append(blks[0].g_db)
    d['nf_c'].append(blks[0].nf_db)
    d['oip3_c'].append(blks[0].oip3_dbm)
    d['oip2_c'].append(blks[0].oip2_dbm)
    d['op1db_c'].append(blks[0].op1db_dbm)

    for k in range(1, len(blks)):
        g_running += blks[k].g_db
        d['gain_c'].append(g_running)
        d['nf_c'].append(calc_nf_db(d['nf_c'][k-1], blks[k].nf_db, d['gain_c'][k-1]))
        d['oip3_c'].append(calc_oip3_dbm(d['oip3_c'][k-1], blks[k].oip3_dbm, blks[k].g_db))
        d['oip2_c'].append(calc_oip2_dbm(d['oip2_c'][k-1], blks[k].oip2_dbm, blks[k].g_db))
        d['op1db_c'].append(calc_op1db_dbm(d['op1db_c'][k-1], blks[k].op1db_dbm, blks[k].g_db))

    return d


def get_cascade_dataframe(blk_df):
    """
    return a Pandas.DataFrame with the entire cascade
    """
    d = {
        'gain_c': [],
        'nf_c': [],
        'oip3_c': [],
        'oip2_c': [],
        'op1db_c': [],
        'iip3_c': [],
        'iip2_c': [],
        'ip1db_c': []
    }
    g_running = blk_df.g_db[0]
    d['gain_c'].append(blk_df.g_db[0])
    d['nf_c'].append(blk_df.nf_db[0])
    d['oip3_c'].append(blk_df.oip3_dbm[0])
    d['oip2_c'].append(blk_df.oip2_dbm[0])
    d['op1db_c'].append(blk_df.op1db_dbm[0])

    d['iip3_c'].append(blk_df.oip3_dbm[0] - blk_df.g_db[0])
    d['iip2_c'].append(blk_df.oip2_dbm[0] - blk_df.g_db[0])
    d['ip1db_c'].append(blk_df.op1db_dbm[0] - blk_df.g_db[0])

    for k in range(1, len(blk_df.index)):
        g_running += blk_df.g_db[k]
        d['gain_c'].append(g_running)
        d['nf_c'].append(calc_nf_db(d['nf_c'][k-1], blk_df.nf_db[k], d['gain_c'][k-1]))
        d['oip3_c'].append(calc_oip3_dbm(d['oip3_c'][k-1], blk_df.oip3_dbm[k], blk_df.g_db[k]))
        d['oip2_c'].append(calc_oip2_dbm(d['oip2_c'][k-1], blk_df.oip2_dbm[k], blk_df.g_db[k]))
        d['op1db_c'].append(calc_op1db_dbm(d['op1db_c'][k-1], blk_df.op1db_dbm[k], blk_df.g_db[k]))

        d['iip3_c'].append(d['oip3_c'][k] - d['gain_c'][k])
        d['iip2_c'].append(d['oip2_c'][k] - d['gain_c'][k])
        d['ip1db_c'].append(d['op1db_c'][k] - d['gain_c'][k])

    for k in d.keys():
        blk_df[k] = pd.Series(d[k], index=blk_df.index)

    # return d
    return blk_df


def calc_op1db_dbm(p1db_up_to_input_dbm, p1db_curr_dbm, g_curr_db):
    """
    return cascaded OIP3 value in dBm after current stage given
    OIP3 up to current stage and OIP3 and gain of current
    stage by itself
    """
    p1db_input = dbm_to_w(p1db_up_to_input_dbm)
    p1db_curr = dbm_to_w(p1db_curr_dbm)
    g_curr = db_to_ratio(g_curr_db)
    p1db_c = 1/(1/(p1db_input*g_curr) + 1/p1db_curr)
    return w_to_dbm(p1db_c)


def calc_oip2_dbm(oip2_up_to_input_dbm, oip2_curr_dbm, g_curr_db):
    """
    return cascaded OIP3 value in dBm after current stage given
    OIP3 up to current stage and OIP3 and gain of current
    stage by itself
    """
    oip2_input = dbm_to_w(oip2_up_to_input_dbm)
    oip2_curr = dbm_to_w(oip2_curr_dbm)
    g_curr = db_to_ratio(g_curr_db)
    # print("oip2_input: {:<0.2f}".format(oip2_input))
    # print("oip2_curr: {:<0.2f}".format(oip2_curr))
    # print("g_curr: {:<0.2f}".format(g_curr))
    ip2_c = 1/(np.sqrt(1/((oip2_input*g_curr)) + np.sqrt(1/oip2_curr))**2)
    return w_to_dbm(ip2_c)


def calc_oip3_dbm(oip3_up_to_input_dbm, oip3_curr_dbm, g_curr_db):
    """
    return cascaded OIP3 value in dBm after current stage given
    OIP3 up to current stage and OIP3 and gain of current
    stage by itself
    """
    oip3_input = dbm_to_w(oip3_up_to_input_dbm)
    oip3_curr = dbm_to_w(oip3_curr_dbm)
    g_curr = db_to_ratio(g_curr_db)
    ip3_c = 1/(1/(oip3_input*g_curr) + 1/oip3_curr)
    return w_to_dbm(ip3_c)


def calc_nf_db(nf_input_db, nf_curr_db, g_up_to_input_db):
    """
    return cascaded noise figure value in dB after current stage given
    noise figure and gain up to current stage and noise figure
    of the current stage by itself
    """
    f_input = db_to_ratio(nf_input_db)
    f_curr = db_to_ratio(nf_curr_db)
    g_up_to_input = db_to_ratio(g_up_to_input_db)
    f = f_input + (f_curr - 1) / g_up_to_input
    return ratio_to_db(f)


def read_cascade_file(fname):
    """

    """
    fo = open(fname, 'r')
    lines = fo.readlines()
    hdr = lines[0].split()
    blks = []
    for line in lines[1:]:
        if line == "\n":
            continue
        d = {}
        for k in range(len(hdr)):
            try:
                d[hdr[k]] = line.split()[k]
            except IndexError:
                break
        blks.append(Block(**d).get_dict())
    return pd.DataFrame(blks)
    # blks.append(Block(**d))
    # return blks


def print_cascade(blks):
    """

    """
    hdr = ['name',
           'gain',
           'nf',
           'oip3',
           'oip2',
           'op1db',
           'gain_c',
           'nf_c',
           'oip3_c',
           'oip2_c',
           'op1db_c']
    print("{:<20} | {:<10}{:<10}{:<10}{:<10}{:<10} | {:<10}{:<10}{:<10}{:<10}{:<10}".format(*hdr))
    d = get_cascade(blks)
    for k in range(len(blks)):
        dat = []
        dat.append(blks[k].name)
        dat.append(blks[k].g_db)
        dat.append(blks[k].nf_db)
        dat.append(blks[k].oip3_dbm)
        dat.append(blks[k].oip2_dbm)
        dat.append(blks[k].op1db_dbm)
        dat.append(d['gain_c'][k])
        dat.append(d['nf_c'][k])
        dat.append(d['oip3_c'][k])  # oip3_c
        dat.append(d['oip2_c'][k])  # oip2_c
        dat.append(d['op1db_c'][k])  # op1db_c
        print("{:<20} | {:<10.1f}{:<10.1f}{:<10.1f}{:<10.1f}{:<10.1f} | {:<10.1f}{:<10.1f}{:<10.1f}{:<10.1f}{:<10.1f}".format(*dat))

    # return d

if __name__=='__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("f", type=str, help="filename")
    args = parser.parse_args()

    blk_df = read_cascade_file(args.f)
    blk_df = get_cascade_dataframe(blk_df)
    print(tabulate(blk_df, headers='keys'))

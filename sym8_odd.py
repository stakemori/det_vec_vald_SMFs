# -*- coding: utf-8; mode: sage -*-
'''
Before running this script,
1. Please install the package degree2 (https://github.com/stakemori/degree2).
2. Please create a directory for saving the genrators of the module.
 The default directory name is "~/sym8_data".
 To change it, please change the 11th line.
'''
import os

sym8_odd_filename_base = os.path.join(os.environ["HOME"], "sym8_data")

from sage.all import fork
from degree2.all import (degree2_modular_forms_ring_level1_gens,
                         Deg2ModularFormQseries)

from degree2.deg2_fourier import SymmetricWeightModularFormElement \
    as SWMFE

from degree2.rankin_cohen_diff import (rankin_cohen_triple_det_sym8,
                                       _rankin_cohen_pair_sym_pol,
                                       rankin_cohen_pair_x5,
                                       vector_valued_rankin_cohen,
                                       rankin_cohen_pair_det2_sym)

from degree2.interpolate import det_deg2

from degree2.utils import mul


sym8_odd_prec = 18

sym8_odd_dct = {"f13": ({4: 1}, {4: 1}, {4: 1}),
                "f15": ({4: 1}, {4: 1}, {6: 1}),
                "g15": ({6: 1}, {4: 1}, {4: 1}),
                "f17": ({4: 1}, {6: 1}, {6: 1}),
                "g17": ({6: 1}, {6: 1}, {4: 1}),
                "f19": ({10: 1}, {4: 1}, {4: 1}),
                "f23": ({4: 1}, {6: 1}, {12: 1})}


def fname_sym8_odd(key, prec):
    return os.path.join(sym8_odd_filename_base, key + "_prec" + str(prec) +
                        ".sobj")


@fork
def save_sym8_odd_form(key, prec=sym8_odd_prec):
    fname = fname_sym8_odd(key, prec)
    if os.path.exists(fname):
        return None
    es4, es6, x10, x12, x35 = degree2_modular_forms_ring_level1_gens(prec)
    if key in sym8_odd_dct:
        wt_dcts = sym8_odd_dct[key]
        gens_dct = {4: es4, 6: es6, 10: x10, 12: x12, 35: x35}
        fs = [mul([gens_dct[k]**i for k, i in d.items()]) for d in wt_dcts]
        f = rankin_cohen_triple_det_sym8(*fs)
    elif key == "h15":
        f10 = rankin_cohen_pair_x5(_rankin_cohen_pair_sym_pol(8, 5, 5), prec)
        f = vector_valued_rankin_cohen(es4, f10)
    elif key == "h17":
        f12 = rankin_cohen_pair_det2_sym(8, es4, es6)
        f = vector_valued_rankin_cohen(es4, f12)
    f.save_as_binary(fname)


def save_all_gens_sym8_odd(prec=sym8_odd_prec):
    keys = ["f13", "f15", "g15", "f17", "g17", "f19", "f23", "h15", "h17"]
    for f in keys:
        save_sym8_odd_form(f, prec=prec)


@fork
def det_sym8_odd(prec=sym8_odd_prec):
    fname = fname_sym8_odd("f187", prec)
    if os.path.exists(fname):
        return None
    keys = ["f13", "f15", "g15", "f17", "g17", "f19", "f23", "h15", "h17"]
    lst = [SWMFE.load_from(fname_sym8_odd(key, prec)) for key in keys]
    f187 = det_deg2([f.forms for f in lst], autom=True,
                    wt=sum([f.wt for f in lst]) + 36)
    f187.save_as_binary(fname)


def check_sym8_odd(prec=sym8_odd_prec):
    det_sym8_odd(prec=sym8_odd_prec)
    f187 = Deg2ModularFormQseries.load_from(fname_sym8_odd("f187", prec))
    es4, es6, _, _, x35 = degree2_modular_forms_ring_level1_gens(prec)
    g187 = (es4**3 - es6**2) * x35**5
    g187 = g187.normalize((12, -1, 14))
    f187 = f187.normalize((12, -1, 14))
    assert f187 == g187

save_all_gens_sym8_odd()
check_sym8_odd()

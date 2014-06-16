# -*- coding: utf-8; mode: sage -*-
'''
Before running this script,
1. Please install the package degree2 (https://github.com/stakemori/degree2).
2. Please create a directory for saving the genrators of the module.
 The default directory name is "~/sym8_data".
 To change it, please change the 11th line.
'''
import os

sym8_even_filename_base = os.path.join(os.environ["HOME"], "sym8_data")

from sage.all import fork
from degree2.all import (rankin_cohen_pair_det2_sym,
                         rankin_cohen_pair_sym,
                         degree2_modular_forms_ring_level1_gens,
                         Deg2ModularFormQseries)

from degree2.deg2_fourier import SymmetricWeightModularFormElement \
    as SWMFE

from degree2.rankin_cohen_diff import (_rankin_cohen_pair_sym_pol,
                                       rankin_cohen_pair_x5)

from degree2.interpolate import det_deg2
from degree2.utils import mul


sym8_even_prec = 15

sym8_even_dct = {"f8": {"weights": ({4: 1}, {4: 1}), "a": 0},
                 "g10": {"weights": ({4: 1}, {6: 1}), "a": 0},
                 "f12": {"weights": ({4: 1}, {6: 1}), "a": 2},
                 "g12": {"weights": ({6: 1}, {6: 1}), "a": 0},
                 "f14": {"weights": ({4: 1}, {10: 1}), "a": 0},
                 "g14": {"weights": ({6: 1}, {4: 2}), "a": 0},
                 "f16": {"weights": ({6: 1}, {10: 1}), "a": 0},
                 "f18": {"weights": ({6: 1}, {12: 1}), "a": 0}}


def fname_sym8_even(key, prec):
    return os.path.join(sym8_even_filename_base,
                        key + "_prec" + str(prec) + ".sobj")


@fork
def save_sym8_even_form(key, prec=sym8_even_prec):
    fname = fname_sym8_even(key, prec)
    if os.path.exists(fname):
        return None

    if key in sym8_even_dct:
        d1, d2 = sym8_even_dct[key]["weights"]
        a = sym8_even_dct[key]["a"]
        es4, es6, x10, x12, x35 = degree2_modular_forms_ring_level1_gens(prec)
        gens_dct = {4: es4, 6: es6, 10: x10, 12: x12, 35: x35}

        f = mul([gens_dct[k]**v for k, v in d1.items()])
        g = mul([gens_dct[k]**v for k, v in d2.items()])

        if a == 0:
            diff_pair_func = rankin_cohen_pair_sym
        elif a == 2:
            diff_pair_func = rankin_cohen_pair_det2_sym

        res = diff_pair_func(8, f, g)
    elif key == "f10":
        res = rankin_cohen_pair_x5(_rankin_cohen_pair_sym_pol(8, 5, 5), prec)
    res.save_as_binary(fname)


def save_all_gens_sym8_even(prec=sym8_even_prec):
    keys = ["f8", "f10", "g10", "f12", "g12", "f14", "g14", "f16", "f18"]
    for f in keys:
        save_sym8_even_form(f, prec=prec)


@fork
def det_sym8_even(prec=sym8_even_prec):
    fname = fname_sym8_even("f150", prec)
    if os.path.exists(fname):
        return None

    keys = ["f8", "f10", "g10", "f12", "g12", "f14", "g14", "f16", "f18"]
    lst = [SWMFE.load_from(fname_sym8_even(n, prec)) for n in keys]
    f150 = det_deg2([f.forms for f in lst], autom=True,
                    wt=sum([f.wt for f in lst]) + 36)
    f150.save_as_binary(fname)


def check_sym8_even(prec=sym8_even_prec):
    _, _, x10, _, x35 = degree2_modular_forms_ring_level1_gens(prec)
    det_sym8_even(prec=prec)
    f150 = Deg2ModularFormQseries.load_from(fname_sym8_even("f150", prec))
    g150 = x10 * x35**4
    assert f150.normalize((9, 5, 13)) == g150


save_all_gens_sym8_even()
check_sym8_even()

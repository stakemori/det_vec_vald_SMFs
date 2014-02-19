# -*- coding: utf-8; mode: sage -*-
from sage.all import matrix

from degree2.utils import naive_det_func
from degree2.all import (
    rankin_cohen_pair_sym, eisenstein_series_degree2,
    x10_with_prec, x12_with_prec, x35_with_prec,
    rankin_cohen_pair_det2_sym,
    Deg2QsrsElement,
    SymmetricWeightModularFormElement,
    Deg2ModularFormQseries,
    CuspFormsDegree2)
from degree2.basic_operation import PrecisionDeg2


def prec_down(f, prec):
    prec = PrecisionDeg2(prec)
    forms_res = []
    for b in f.forms:
        d = b._to_format_dct()
        d["prec"] = prec._to_format_dct()
        b = Deg2QsrsElement._from_dict_to_object(d)
        forms_res.append(b)
    return SymmetricWeightModularFormElement(forms_res, f.wt, prec)

# For a scalar valued Siegel modular form f of degree2, level 1,
# weight 115, if a((n, r, m), f) = 0 for (n, r, m) in tpls115 below,
# then f = 0.
tpls115 = [(3, 1, 7),
           (6, 4, 8),
           (2, 1, 9),
           (5, 3, 9),
           (4, 1, 6),
           (3, 1, 4),
           (6, 5, 7),
           (4, 2, 6),
           (4, 2, 8),
           (6, 4, 9),
           (7, 1, 11),
           (3, 2, 6),
           (2, 1, 8),
           (7, 4, 11),
           (5, 3, 6),
           (3, 2, 5),
           (5, 1, 6),
           (4, 2, 10),
           (6, 1, 10),
           (6, 1, 11),
           (5, 4, 9),
           (4, 3, 7),
           (4, 3, 6),
           (5, 3, 8),
           (3, 2, 4),
           (9, 8, 10),
           (3, 1, 5),
           (4, 3, 5),
           (2, 1, 6),
           (3, 2, 8),
           (10, 2, 11),
           (5, 4, 7),
           (3, 2, 7),
           (7, 6, 8),
           (4, 3, 10),
           (8, 7, 9),
           (10, 9, 11),
           (8, 7, 10),
           (6, 3, 7),
           (6, 4, 10),
           (6, 3, 8),
           (2, 1, 10),
           (4, 1, 5),
           (5, 2, 8),
           (7, 5, 11),
           (3, 1, 8),
           (4, 2, 5),
           (6, 2, 7),
           (8, 4, 11),
           (9, 7, 11),
           (3, 1, 6),
           (4, 3, 8),
           (7, 6, 9),
           (6, 5, 8),
           (5, 1, 10),
           (2, 1, 7),
           (4, 1, 10),
           (4, 3, 9),
           (5, 3, 10),
           (7, 5, 9),
           (7, 4, 9),
           (8, 7, 11),
           (5, 2, 11),
           (5, 4, 8),
           (2, 1, 5),
           (3, 2, 10),
           (6, 1, 9),
           (9, 4, 11),
           (3, 1, 10),
           (8, 2, 9),
           (6, 3, 10),
           (7, 1, 8),
           (8, 4, 10),
           (5, 1, 9),
           (2, 1, 3),
           (6, 5, 10),
           (5, 3, 11),
           (6, 4, 7),
           (5, 4, 6),
           (8, 1, 9),
           (8, 6, 9),
           (8, 2, 10),
           (8, 3, 11),
           (4, 1, 7),
           (10, 8, 11),
           (6, 2, 11),
           (4, 2, 7),
           (6, 5, 9),
           (10, 3, 11),
           (10, 7, 11),
           (5, 3, 7),
           (7, 3, 8),
           (5, 2, 7),
           (5, 4, 10),
           (9, 4, 10),
           (4, 2, 11),
           (9, 5, 10),
           (6, 2, 10),
           (7, 2, 10),
           (10, 1, 11),
           (7, 3, 11)]

S115 = CuspFormsDegree2(115)
assert matrix([[h[t] for h in S115.basis()]
               for t in tpls115]).rank() == S115.dimension()


# weight 115 Siegel modualr form is equal to 0 if
# the Fourier coefficients are 0 in this precision.
prec_wt115 = PrecisionDeg2(tpls115)

es4 = eisenstein_series_degree2(4, prec_wt115)
es6 = eisenstein_series_degree2(6, prec_wt115)
x10 = x10_with_prec(prec_wt115)
x12 = x12_with_prec(prec_wt115)

# Siegel Eisenstein series of degree 2 with precision 22
# and weight 4 and 6 respectively.
es4_prec22 = eisenstein_series_degree2(4, prec=22)
es6_prec22 = eisenstein_series_degree2(6, prec=22)

F12_prec22 = rankin_cohen_pair_det2_sym(6, es4_prec22, es6_prec22)

g12_fc_dct = {}
for t in prec_wt115:
    g12_fc_dct[t] = F12_prec22.hecke_operator(2, t)

g12_forms = []
for i in range(7):
    dct = {t: v.vec[i] for t, v in g12_fc_dct.iteritems()}
    f1 = Deg2QsrsElement(dct, prec_wt115, is_cuspidal=True)
    g12_forms.append(f1)

F12 = prec_down(F12_prec22, prec_wt115)
G12 = SymmetricWeightModularFormElement(g12_forms, 12, prec_wt115)
H12 = rankin_cohen_pair_sym(6, es4, es4**2)
F10 = rankin_cohen_pair_sym(6, es4, es6)
F14 = rankin_cohen_pair_sym(6, es4, x10)
F16 = rankin_cohen_pair_sym(6, es4, x12)
F18 = rankin_cohen_pair_sym(6, es6, x12)

lst = [F12, G12, H12, F10, F14, F16, F18]

det_form = naive_det_func(7)([a.forms for a in lst])
det_form = Deg2ModularFormQseries(
    115, det_form.fc_dct, det_form.prec, is_cuspidal=True)

x35 = x35_with_prec(prec_wt115)
wt115form = x35**3 * es4 * es6

assert wt115form.normalize((7, -1, 8)) == det_form.normalize((7, -1, 8))

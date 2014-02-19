# -*- coding: utf-8; mode: sage -*-
from sage.all import QQ

from degree2.utils import naive_det_func
from degree2.all import (
    rankin_cohen_pair_sym, eisenstein_series_degree2,
    x10_with_prec, x12_with_prec, x35_with_prec,
    rankin_cohen_pair_det2_sym,
    rankin_cohen_triple_det_sym4)

det_5 = naive_det_func(5)
prec_ev = 7

phi4 = eisenstein_series_degree2(4, prec_ev)
phi6 = eisenstein_series_degree2(6, prec_ev)
x10 = x10_with_prec(prec_ev)
x12 = x12_with_prec(prec_ev)
x35 = x35_with_prec(prec_ev)

f1 = rankin_cohen_pair_sym(4, phi4, phi4)
f2 = rankin_cohen_pair_sym(4, phi4, phi6)
f3 = rankin_cohen_pair_det2_sym(4, phi4, phi6)
f4 = rankin_cohen_pair_sym(4, phi4, x10)
f5 = rankin_cohen_pair_sym(4, phi6, x10)

x70 = det_5([f.forms for f in [f1, f2, f3, f4, f5]])
y70 = x70 * QQ(-19945421021123916595200000)**(-1)
assert y70 == x35**2


prec_od = 10

e4 = eisenstein_series_degree2(4, prec_od)
e6 = eisenstein_series_degree2(6, prec_od)
x10 = x10_with_prec(prec_od)
x12 = x12_with_prec(prec_od)
x35 = x35_with_prec(prec_od)

g0 = rankin_cohen_triple_det_sym4(e4, e4, e6)
g1 = rankin_cohen_triple_det_sym4(e4, e6, e6)
g2 = rankin_cohen_triple_det_sym4(e4, e4, x10)
g3 = rankin_cohen_triple_det_sym4(e4, e4, x12)
g4 = rankin_cohen_triple_det_sym4(e4, e6, x12)

x105 = det_5([f.forms for f in [g0, g1, g2, g3, g4]])
y105 = x105 * QQ(75141695851970512394649600000000)**(-1)
assert y105 == x35**3

# -*- coding: utf-8; mode: sage -*-

from degree2.all import (eisenstein_series_degree2, x10_with_prec,
                         x12_with_prec, x35_with_prec,
                         rankin_cohen_pair_sym,
                         rankin_cohen_pair_det2_sym, Deg2ModularFormQseries,
                         Deg2QsrsElement)


from sage.all import PolynomialRing, QQ, matrix

from degree2.utils import naive_det

from degree2.rankin_cohen_diff import (
    rankin_cohen_triple_x5,
    _rankin_cohen_bracket_func)

from degree2.basic_operation import PrecisionDeg2

from degree2.deg2_fourier import SymmetricWeightGenericElement\
    as SWGElt
from degree2.deg2_fourier import SymmetricWeightModularFormElement \
    as SWMFE

from degree2.deg2_fourier import _common_base_ring

prec_od = 15

def _odd_sym6_pol_dct():
    R = PolynomialRing(QQ, names="r,s,t")
    r, s, t = R.gens()
    dnm = {15: QQ(160),
           17: QQ(192),
           19: QQ(1920),
           21: QQ(2880),
           23: QQ(16)}
    nm = {15: (5, -14, 7),
          17: (4, -8, 3),
          19: (22, -24, 5),
          21: (22, -24, 5),
          23: (13, -14, 3)}
    res = {}
    for k in dnm.keys():
        a, b, c = nm[k]
        res[k] = (a*r**2 + b * r * t + c*t**2) * QQ(1)/QQ(dnm[k])
    return res

def cross_prod(v1, v2):
    a, b, c = v1
    ad, bd, cd = v2

    return (2 * (a * bd - b * ad),
            a * cd - c * ad,
            2 * (b * cd - c * bd))


def m_operator(k1, k2, k3, rnames=None, unames=None):
    '''The operator M_k by van Dorp.
    '''
    if rnames is None:
        rnames = "r11, r12, r22, s11, s12, s22, t11, t12, t22"
    if unames is None:
        unames = "u1, u2"
    R = PolynomialRing(QQ, names=rnames)
    S = PolynomialRing(R, names=unames)
    r11, r12, r22, s11, s12, s22, t11, t12, t22 = R.gens()
    rs = (r11, r12, r22)
    ss = (s11, s12, s22)
    ts = (t11, t12, t22)
    u1, u2 = S.gens()

    def bracket_op(rs):
        r1, r2, r3 = rs
        return r1 * u1**2 + 2 * r2*u1*u2 + r3 * u2**2

    def x_op_val(f):
        r, s, t = f.parent().gens()
        return f.subs({r: bracket_op(rs),
                       s: bracket_op(ss),
                       t: bracket_op(ts)})

    def m_op_val(f):
        r, s, t = f.parent().gens()
        x_val = x_op_val(f)
        xs = [k * x_val for k in [k3, k2, k1]]
        brxs = [bracket_op(a) * x_op_val(f.derivative(b))
                for a, b in zip([ts, ss, rs], [t, s, r])]
        brcks = [bracket_op(cross_prod(a, b))
                 for a, b in zip([rs, ts, ss], [ss, rs, ts])]
        return sum([a * (b + c) for a, b, c in zip(brcks, xs, brxs)])

    return m_op_val

def diff_u(vec_val):
    forms = [i * f for f, i in zip(vec_val.forms,
                                   reversed(range(1, vec_val.sym_wt + 1)))]
    return SWGElt(forms, vec_val.prec, vec_val.base_ring)

def diff_v(vec_val):
    forms = [i * f for f, i in zip(vec_val.forms[1:],
                                   range(1, vec_val.sym_wt + 1))]
    return SWGElt(forms, vec_val.prec, vec_val.base_ring)

def diff_d(vec_val):
    return [diff_u(diff_u(vec_val)),
            diff_u(diff_v(vec_val)),
            diff_v(diff_v(vec_val))]


def differential_monom(vec_val, a, b, c):
    forms = [f._differential_operator_monomial(a, b, c)
             for f in vec_val.forms]
    return SWGElt(forms, vec_val.prec, vec_val.base_ring)


def cross_prod_diff(vec_vals):
    f1, f2, f3 = vec_vals

    def d1(f):
        return differential_monom(f, 1, 0, 0)

    def d2(f):
        return differential_monom(f, 0, 1, 0) * QQ(2)**(-1)

    def d3(f):
        return differential_monom(f, 0, 0, 1)

    return [2 * (d1(f2) - d2(f1)),
            d1(f3) - d3(f1),
            2 * (d2(f3) - d3(f2))]

def bracket_vec_val(vecs):
    if isinstance(vecs[0], SWGElt):
        v1, v2, v3 = [a.forms for a in vecs]
    else:
        v1, v2, v3 = vecs
    j = len(v1) - 1
    def _names(s):
        return ", ".join([s + str(i) for i in range(j + 1)])

    R = PolynomialRing(QQ, names=", ".join([_names(s) for s in
                                            ["x", "y", "z"]]))
    gens_x = R.gens()[: j + 1]
    gens_y = R.gens()[j + 1 : 2 * (j + 1)]
    gens_z = R.gens()[2 * (j + 1):]
    S = PolynomialRing(R, names="u, v")
    u, v = S.gens()

    def _pol(gens):
        return sum([a * u**(j - i) * v**i
                    for i, a in zip(range(j + 1), gens)])

    f_x, f_y, f_z = [_pol(gens) for gens in [gens_x, gens_y, gens_z]]
    A = matrix([[f_x, f_y],
                [f_y, f_z]])
    vec = matrix([u, v]).transpose()
    g = (vec.transpose() * A * vec)[0][0]
    pol_dc = {(i, j + 2 - i): g[(i, j + 2 - i)] for i in range(j + 3)}

    def pol_to_val(f):
        dct = {}
        def _dct(gens, v):
            return {a: b for a, b in zip(gens, v)}

        dct.update(_dct(gens_x, v1))
        dct.update(_dct(gens_y, v2))
        dct.update(_dct(gens_z, v3))
        return f.subs(dct)

    res_dc = {k: pol_to_val(v) for k, v in pol_dc.iteritems()}
    return [res_dc[(j + 2 - i, i)] for i in range(j + 3)]



def vector_valued_rankin_cohen(f, vec_val):
    '''
    Rankin-Cohen type bracket defined by van Dorp.
    '''
    sym_wt = vec_val.sym_wt
    prec = f.prec
    base_ring = _common_base_ring(f.base_ring, vec_val.base_ring)
    diff_tau = (f.differentiate_wrt_tau(),
                f.differentiate_wrt_z() * QQ(2)**(-1),
                f.differentiate_wrt_w())
    crs_prd1 = cross_prod(diff_tau, diff_d(vec_val))
    forms1 = bracket_vec_val(crs_prd1)
    res1 = (vec_val.wt + sym_wt//2 - 1) * SWGElt(forms1, prec,
                                           base_ring=base_ring)

    forms2 = bracket_vec_val(cross_prod_diff(diff_d(vec_val)))
    res2 = f.wt * f * SWGElt(forms2, prec, base_ring=base_ring)

    res = SWMFE((res1 - res2).forms, f.wt + vec_val.wt + 1,
                prec, base_ring=base_ring)
    return res

def _sym6_od_Q(k, t):
    return m_operator(*t)(_odd_sym6_pol_dct()[k])

def f15_sym6(prec):
    '''One of the generator with the determinant weight 15 given by
    van Dorp.
    '''
    Q = _sym6_od_Q(15, (5, 5, 4))
    phi4 = eisenstein_series_degree2(4, prec + 1)
    return rankin_cohen_triple_x5(Q, phi4, prec)

def f17_sym6(prec):
    '''One of the generator with the determinant weight 17 given by
    van Dorp.
    '''
    Q = _sym6_od_Q(17, (5, 5, 6))
    phi6 = eisenstein_series_degree2(6, prec + 1)
    return rankin_cohen_triple_x5(Q, phi6, prec)

def f19_sym6(prec):
    '''One of the generator with the determinant weight 19 given by
    van Dorp.
    '''
    Q = _sym6_od_Q(19, (4, 4, 10))
    phi4 = eisenstein_series_degree2(4, prec)
    x10 = x10_with_prec(prec)
    forms = _rankin_cohen_bracket_func(Q)([phi4, phi4, x10])
    return SWMFE(forms, 19, prec)

def f21_sym6(prec):
    '''One of the generator with the determinant weight 21 given by
    van Dorp.
    '''
    Q = _sym6_od_Q(21, (4, 6, 10))
    phi4 = eisenstein_series_degree2(4, prec)
    phi6 = eisenstein_series_degree2(6, prec)
    x10 = x10_with_prec(prec)
    forms = _rankin_cohen_bracket_func(Q)([phi4, phi6, x10])
    return SWMFE(forms, 21, prec)

def f23_sym6(prec):
    '''One of the generator with the determinant weight 23 given by
    van Dorp.
    '''
    Q = _sym6_od_Q(23, (5, 5, 12))
    x12 = x12_with_prec(prec + 1)
    return rankin_cohen_triple_x5(Q, x12, prec)

f15 = f15_sym6(prec_od)
f17 = f17_sym6(prec_od)
f19 = f19_sym6(prec_od)
f21 = f21_sym6(prec_od)
f23 = f23_sym6(prec_od)

# The rank of A_{Sym(6)}^{1}(Gamma_2) is equal to 7.
# Therefore we need other two vector valued Siegel modular forms.

# We construct two modular forms f17_1 and f17_2 by
# Siegel modualr forms of weight det^12 Sym(6) and
# Rankin-Cohen type differential operator defined by Dorp
# (i.e. function vector_valued_rankin_cohen).

es4 = eisenstein_series_degree2(4, prec_od)
es6 = eisenstein_series_degree2(6, prec_od)
x35 = x35_with_prec(prec_od)
f12 = rankin_cohen_pair_det2_sym(6, es4, es6)
h12 = rankin_cohen_pair_sym(6, es4, es4**2)

f17_1 = vector_valued_rankin_cohen(es4, f12)

# f12_prect2 is {phi4, phi6}_{det^2 Sym(6)}.
# The prec of f12_prect2 is larger than other forms.
f12_prect2 = rankin_cohen_pair_det2_sym(
    6,
    eisenstein_series_degree2(4, prec_od*2),
    eisenstein_series_degree2(6, prec_od*2))

def _g12():
    g12_fc_dct = {}
    for t in PrecisionDeg2(prec_od):
        g12_fc_dct[t] = f12_prect2.hecke_operator(2, t)
    g12_forms = []
    for i in range(7):
        dct = {t: v.vec[i] for t, v in g12_fc_dct.iteritems()}
        f1 = Deg2QsrsElement(dct, prec_od, is_cuspidal=True)
        g12_forms.append(f1)
    return SWMFE(g12_forms, 12, prec_od)

# g12 is the image of f12_prect2 by the Hecke operator T(2).
g12 = _g12()

f17_2 = vector_valued_rankin_cohen(es4, g12)

# The determinant of [f15, f17, f17_1, f17_2, f19, f21, f23]
f150 = Deg2ModularFormQseries(
    150,
    naive_det([f15, f17, f17_1, f17_2, f19, f21, f23]),
    prec_od,
    is_cuspidal=True)

# Check the determinant is equal to phi4 phi6 x35^4 up to a constant.
assert f150.normalize((8, -4, 12)) == es4 * es6 * x35**4

from degree2.vector_valued_impl.sym10.odd_structure import calculator
from degree2.interpolate import det_deg2
from degree2.scalar_valued_smfs import eisenstein_series_degree2, x35_with_prec

def check_det_with_prec(prec):
    es4 = eisenstein_series_degree2(4, prec)
    es6 = eisenstein_series_degree2(6, prec)
    x35 = x35_with_prec(prec)
    f = (es4**3 - es6**2) * x35**6
    t = f._none_zero_tpl()
    d = calculator.forms_dict(prec)
    # Constructions of the first 11 generators.
    cs = list(sorted(calculator._const_vecs, key=lambda x: x.weight()))[:11]
    wt = sum(c.weight() for c in cs) + (10 * 11)//2
    mat = [d[c].forms for c in cs]
    det = det_deg2(mat, wt=wt)
    assert f[t] * det == det[t] * f

check_det_with_prec(22)

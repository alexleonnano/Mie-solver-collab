import numpy as np

def wvl_interpolation(wvl_range, ref_lam, ref_n, ref_k):
    # Interpolate refractive index
    n = np.interp(wvl_range, ref_lam, ref_n)
    k = np.interp(wvl_range, ref_lam, ref_k)

    return n, k
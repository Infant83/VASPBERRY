"""Occupation-aware contractions of rectangular standard optical connections.

No square matrix is reconstructed. For i<j, K_ij=-2 Im conj(C_ji,a) C_ji,b;
its coefficient is w_i-w_j, with w_n=1[n in S] for a geometric trace or
w_n=1[n in S]*f(E_n,mu,T) for conductivity. S selects target states, not
the virtual-state cutoff: pair geometry retains every source NBANDS state.
C already includes its energy denominator, so no second division applies.
Producer-erased and absent pairs can cancel only when their physical
response weights are equal.
"""
from __future__ import annotations

from types import SimpleNamespace
import numpy as np

from berry_data import require
from band_selection import parse_bands
from vaspberry_transport import KB_EV_PER_K, fermi_dirac

PRODUCER_THRESHOLD_EV = .002
DENOMINATOR_THRESHOLD_EV = float(np.float32(1e-10))



def selected_mask(bands, nb):
    require(isinstance(bands, (list, tuple, np.ndarray)) and len(bands) > 0,
            'nonempty selected bands required')
    require(all(type(v) in (int, np.int64, np.int32) and 1 <= v <= nb for v in bands)
            and len(set(bands)) == len(bands), 'selected bands must be distinct source IDs within NBANDS')
    return np.isin(np.arange(1, nb+1), bands)


def optical_pairs(connection, energies):
    """Prepare pair geometry and validity; validate only required pairs later."""
    c = np.asarray(connection); e = np.asarray(energies, dtype=float)
    require(c.ndim == 4 and c.shape[1] == 3 and e.shape == (c.shape[0], c.shape[2])
            and 0 < c.shape[3] <= c.shape[2], 'rectangular optical connection/energy dimensions disagree')
    require(np.isfinite(e).all() and np.all(np.diff(e, axis=1) >= -1e-8), 'finite energy-ordered source bands required')
    nk, _, nb, nd = c.shape
    i, j = np.triu_indices(nb, 1)
    stored = (i < nd) | (j < nd)
    both = (i < nd) & (j < nd)
    elements = np.zeros((nk, 3, len(i)), dtype=np.complex128)
    direct = i < nd
    elements[:, :, direct] = c[:, :, j[direct], i[direct]].astype(np.complex128)
    reverse = stored & ~direct
    elements[:, :, reverse] = c[:, :, i[reverse], j[reverse]].conj().astype(np.complex128)
    finite = np.isfinite(elements).all(axis=1)
    hermitian = np.ones((nk, len(i)), dtype=bool)
    residual_ratio = np.zeros((nk, len(i)))
    gaps = abs(e[:, i]-e[:, j])
    if both.any():
        back = c[:, :, i[both], j[both]].astype(np.complex128).conj()
        forward = elements[:, :, both]
        residual = abs(forward-back)*gaps[:, None, both]
        limit = 1e-10+1e-8*np.maximum(abs(forward), abs(back))*gaps[:, None, both]
        finite[:, both] &= np.isfinite(back).all(axis=1)
        hermitian[:, both] = np.all(residual <= limit, axis=1)
        residual_ratio[:, both] = np.max(residual/limit, axis=1)
    cluster = np.concatenate((np.zeros((nk, 1), dtype=int),
                              np.cumsum(np.diff(e, axis=1) > PRODUCER_THRESHOLD_EV, axis=1)), axis=1)
    erased = (cluster[:, i] == cluster[:, j]) | (gaps <= DENOMINATOR_THRESHOLD_EV)
    omega = np.stack([-2*np.imag(elements[:, a].conj()*elements[:, b])
                      for a, b in ((1, 2), (2, 0), (0, 1))], axis=-1)
    return SimpleNamespace(energies=e, i=i, j=j, omega=omega, stored=stored,
        finite=finite, hermitian=hermitian, residual_ratio=residual_ratio,
        erased=erased, gaps=gaps, clusters=cluster, nbands=nb, ndbands=nd)


def contract_pair_weights(pairs, difference, required=None):
    """Contract one response's pair differences with strict required-pair checks."""
    d = np.asarray(difference, dtype=float)
    require(d.shape == pairs.gaps.shape and np.isfinite(d).all(), 'finite (k,pair) weight differences required')
    needed = d != 0 if required is None else np.asarray(required, dtype=bool)
    require(needed.shape == d.shape and np.all(needed | (d == 0)), 'required-pair mask omits a nonzero weight')
    for invalid, label in ((~pairs.stored[None, :], 'missing rectangular matrix pair'),
                           (pairs.erased, 'producer-erased pair (full-source 0.002 eV cluster/1e-10 eV guard)'),
                           (~pairs.finite, 'nonfinite required optical pair'),
                           (~pairs.hermitian, 'required optical pair violates Hermitian consistency')):
        bad = needed & invalid
        if bad.any():
            k, p = np.argwhere(bad)[0]
            raise ValueError(f'{label}: k={k+1}, bands={pairs.i[p]+1},{pairs.j[p]+1}; unequal response weights')
    # Exclude canceled invalid elements before multiplication, never NaN*0.
    out = np.einsum('kp,kpc->kc', d, np.where(needed[:, :, None], pairs.omega, 0.))
    require(np.isfinite(out).all(), 'nonfinite weighted optical curvature')
    return out


def geometric_trace(pairs, bands):
    """Unit-weight selected trace; complete internal clusters cancel exactly."""
    selected = selected_mask(bands, pairs.nbands).astype(float)
    diff = np.broadcast_to(selected[pairs.i]-selected[pairs.j], pairs.gaps.shape)
    return contract_pair_weights(pairs, diff)


def logistic_difference(x, y):
    """Stable f(x)-f(y), including differences in exponentially small tails."""
    x, y = np.broadcast_arrays(np.asarray(x, float), np.asarray(y, float))
    result = np.empty(x.shape)
    positive = (x >= 0) & (y >= 0); negative = (x < 0) & (y < 0)
    for mask, sign in ((positive, -1.), (negative, 1.)):
        a, b = np.exp(sign*x[mask]), np.exp(sign*y[mask])
        result[mask] = (a-b if sign < 0 else b-a)/((1+a)*(1+b))
    other = ~(positive | negative)
    if other.any():
        def f(z):
            t = np.exp(-abs(z)); return np.where(z >= 0, t/(1+t), 1/(1+t))
        result[other] = f(x[other])-f(y[other])
    return result


def occupation_response(pairs, bands, mu, temperature, mu_reference):
    """Weighted curvature and Delta; exact equality controls missing data."""
    selected = selected_mask(bands, pairs.nbands)
    e = pairs.energies; i, j = pairs.i, pairs.j
    occ = fermi_dirac(e, mu, temperature)*selected
    ref = fermi_dirac(e, mu_reference, temperature)*selected
    if temperature == 0:
        delta_occ = occ-ref
        diff = occ[:, i]-occ[:, j]; refdiff = ref[:, i]-ref[:, j]
        needed = (diff != 0) | (refdiff != 0)
    else:
        # Even underflowed or rounded Fermi tails are physically unequal unless
        # both bands are outside S or both selected energies are exactly equal.
        needed = np.broadcast_to(selected[i] | selected[j], pairs.gaps.shape).copy()
        needed &= ~(selected[i][None, :] & selected[j][None, :] & (e[:, i] == e[:, j]))
        x = (e-mu)/(KB_EV_PER_K*temperature); y = (e-mu_reference)/(KB_EV_PER_K*temperature)
        delta_occ = logistic_difference(x, y)*selected
        diff = occ[:, i]-occ[:, j]
        both = selected[i] & selected[j]
        diff[:, both] = logistic_difference(x[:, i[both]], x[:, j[both]])
    delta_diff = delta_occ[:, i]-delta_occ[:, j]
    return (contract_pair_weights(pairs, diff, needed),
            contract_pair_weights(pairs, delta_diff, needed), occ.sum(axis=1), delta_occ.sum(axis=1), needed)

"""Cleaned reference from the supplied Python algorithm.

Duplicate imports / superseded pure-bisection inverse are removed. Numerical
formulas, default direction and fit decisions are preserved. leggauss tables
are memoized; this does not change the quadrature formula.
"""
from functools import lru_cache
from bisect import bisect_left
from math import atan2, cos, sin, pi, sqrt, ceil
import numpy as np


def build_trajectory_geometry(A, O, R, branch=1, continuous=False):
    Ax, Ay = A
    Ox, Oy = O
    ux, uy = Ax - Ox, Ay - Oy
    D = sqrt(ux * ux + uy * uy)
    if R <= 0 or D <= R + 1e-12:
        raise ValueError("A must be outside a positive-radius circle")
    ex, ey = ux / D, uy / D
    nx, ny = -ey, ex
    ct = R / D
    st = sqrt(max(0.0, 1.0 - ct * ct))
    sign = 1 if branch >= 0 else -1
    px = R * (ct * ex + sign * st * nx)
    py = R * (ct * ey + sign * st * ny)
    Bx, By = Ox + px, Oy + py
    ABx, ABy = Bx - Ax, By - Ay
    s_ab = sqrt(ABx * ABx + ABy * ABy)
    direction = 1.0 if (-py * ABx + px * ABy) < 0 else -1.0
    if continuous:
        direction = -direction
    Cx, Cy = Bx + ABx, By + ABy
    return dict(A=np.array(A), O=np.array(O), R=R, B=np.array([Bx, By]),
                C=np.array([Cx, Cy]), direction=direction, theta0=atan2(py, px),
                s_ab=s_ab, s_circle=2 * pi * R, s_total=2 * s_ab + 2 * pi * R,
                AB=np.array([ABx, ABy]), AC=np.array([Cx - Ax, Cy - Ay]))


def generate_trajectory_points(A, O, R, d, branch=1, continuous=False):
    if not np.isfinite(d) or d <= 0:
        raise ValueError("d must be positive and finite")
    geo = build_trajectory_geometry(A, O, R, branch, continuous)
    sab, sc, total = geo['s_ab'], geo['s_circle'], geo['s_total']
    n = max(1, int(ceil(total / d - 1e-12)))
    if n >= 1000000:
        raise ValueError("control-point limit exceeded")
    step = total / n
    points = []
    for i in range(n + 1):
        s = min(i * step, total)
        if s <= sab:
            pt = geo['A'] + (s / sab) * geo['AB']
        elif s <= sab + sc:
            angle = geo['theta0'] + geo['direction'] * (s - sab) / R
            pt = geo['O'] + R * np.array([cos(angle), sin(angle)])
        else:
            pt = geo['B'] + ((s - sab - sc) / sab) * geo['AB']
        points.append(pt)
    return np.array(points), step, geo


def project_to_line(P, S1, S2):
    v = S2 - S1
    denom = np.dot(v, v)
    if denom < 1e-15:
        return S1.copy(), np.linalg.norm(P - S1)
    Q = S1 + np.dot(P - S1, v) / denom * v
    return Q, np.linalg.norm(P - Q)


def project_to_circle(P, O, R):
    v = P - O
    dist = np.linalg.norm(v)
    if dist < 1e-15:
        return O + R * np.array([1.0, 0.0]), R
    return O + R * v / dist, abs(dist - R)


def project_to_theoretical(P, geo):
    q1, d1 = project_to_line(P, geo['A'], geo['C'])
    q2, d2 = project_to_circle(P, geo['O'], geo['R'])
    return (q1, d1) if d1 <= d2 else (q2, d2)


def find_span(n, p, u, U):
    if u >= U[n + 1]:
        return n
    if u <= U[p]:
        return p
    low, high = p, n + 1
    mid = (low + high) // 2
    while u < U[mid] or u >= U[mid + 1]:
        if u < U[mid]:
            high = mid
        else:
            low = mid
        mid = (low + high) // 2
    return mid


def de_boor(u, P, p, U):
    k = find_span(len(P) - 1, p, u, U)
    work = [P[j].copy() for j in range(k - p, k + 1)]
    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            idx = k - p + j
            denom = U[idx + p - r + 1] - U[idx]
            alpha = (u - U[idx]) / denom if denom > 1e-12 else 0.0
            work[j] = (1 - alpha) * work[j - 1] + alpha * work[j]
    return work[p]


def evaluate_many(us, P, p, U):
    """Vectorized de Boor for dense independent trapezoid checks / plotting."""
    us = np.asarray(us)
    k = np.clip(np.searchsorted(U, us, side='right') - 1, p, len(P) - 1)
    work = np.array([P[k - p + j] for j in range(p + 1)])
    for r in range(1, p + 1):
        for j in range(p, r - 1, -1):
            idx = k - p + j
            denom = U[idx + p - r + 1] - U[idx]
            alpha = np.divide(us - U[idx], denom, out=np.zeros_like(us), where=denom > 1e-12)
            work[j] = (1 - alpha[:, None]) * work[j - 1] + alpha[:, None] * work[j]
    return work[p]


def generate_knot_vector(P, p):
    N = len(P) - 1
    dist = np.linalg.norm(np.diff(P, axis=0), axis=1)
    total = np.sum(dist)
    t = np.linspace(0, 1, N + 1) if total < 1e-12 else np.concatenate(([0], np.cumsum(dist) / total))
    U = np.zeros(N + p + 2)
    U[-(p + 1):] = 1
    for j in range(p + 1, N + 1):
        U[j] = np.sum(t[j - p:j]) / p
    return np.maximum.accumulate(U)


def evaluate_error_for_d(d, A, O, R, branch=1, p=4, sample_num=1000, continuous=False):
    P, step, geo = generate_trajectory_points(A, O, R, d, branch, continuous)
    if len(P) < p + 1:
        return None
    U = generate_knot_vector(P, p)
    curve = np.array([de_boor(u, P, p, U) for u in np.linspace(0, 1, sample_num)])
    errs = [project_to_theoretical(pt, geo)[1] for pt in curve]
    return max(errs), P, U, geo, curve


def fit_bspline_by_bisection(A, O, R, err_min, err_max, branch=1, d_init=.5,
                              max_iter=60, tol_d=1e-5, p=4, sample_num=1000,
                              continuous=False):
    if not 0 < err_min <= err_max:
        raise ValueError("invalid error interval")
    def ev(d):
        res = evaluate_error_for_d(d, A, O, R, branch, p, sample_num, continuous)
        if res is None:
            raise ValueError("insufficient control points")
        return res
    res = ev(d_init)
    if err_min <= res[0] <= err_max:
        return (d_init, *res)
    if res[0] > err_max:
        hi, lo = d_init, d_init / 2
        for _ in range(64):
            if ev(lo)[0] <= err_max:
                break
            hi, lo = lo, lo / 2
        else:
            raise ValueError("failed to bracket")
    else:
        lo, hi = d_init, d_init * 2
        for _ in range(64):
            if ev(hi)[0] >= err_min:
                break
            lo, hi = hi, hi * 2
        else:
            raise ValueError("failed to bracket")
    for _ in range(max_iter):
        mid = (lo + hi) / 2
        res = ev(mid)
        if err_min <= res[0] <= err_max:
            return (mid, *res)
        if res[0] > err_max:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol_d:
            break
    raise ValueError("discrete search did not meet the requested interval")


def compute_derivative_curve(P, p, U):
    if p == 0:
        return None, None
    Q = np.zeros((len(P) - 1, 2))
    for i in range(len(Q)):
        denom = U[i + p + 1] - U[i + 1]
        if denom > 1e-12:
            Q[i] = p * (P[i + 1] - P[i]) / denom
    return Q, U[1:-1]


def bspline_derivative(u, Q, p_minus_1, U_prime):
    return de_boor(u, Q, p_minus_1, U_prime)


@lru_cache(maxsize=16)
def _rule(n):
    return np.polynomial.legendre.leggauss(n)


def gauss_legendre_integral(f, a, b, n):
    nodes, weights = _rule(n)
    half, center = .5 * (b - a), .5 * (b + a)
    return half * sum(weights[i] * f(half * nodes[i] + center) for i in range(n))


class BSplineArcLength:
    """Original query-dependent prefix cache, retained as a comparison baseline."""
    def __init__(self, P, p, U, err=1e-8, n_start=4, n_max=256):
        self.P, self.p, self.U = P, p, U
        self.err, self.n_start, self.n_max = err, n_start, n_max
        self.Q, self.U_prime = compute_derivative_curve(P, p, U)
        self.unique_knots = [U[0]]
        for u in U[1:]:
            if u > self.unique_knots[-1] + 1e-14:
                self.unique_knots.append(u)
        self.cache = [(U[0], 0.)]

    def speed(self, u):
        return float(np.linalg.norm(de_boor(u, self.Q, self.p - 1, self.U_prime)))

    def _integrate(self, a, b, err):
        if b <= a + 1e-15:
            return 0.
        bps = [a] + [u for u in self.unique_knots if a + 1e-14 < u < b - 1e-14] + [b]
        total = 0.
        for aa, bb in zip(bps[:-1], bps[1:]):
            n = self.n_start
            prev = gauss_legendre_integral(self.speed, aa, bb, n)
            while True:
                n2 = 2 * n
                curr = gauss_legendre_integral(self.speed, aa, bb, n2)
                if abs(curr - prev) < err / (len(bps) - 1) or n2 >= self.n_max:
                    total += curr
                    break
                prev, n = curr, n2
        return total

    def cumulative(self, u):
        if u <= self.cache[0][0] + 1e-15:
            return 0.
        idx = bisect_left([pt[0] for pt in self.cache], u - 1e-12)
        if idx < len(self.cache) and abs(self.cache[idx][0] - u) < 1e-12:
            return self.cache[idx][1]
        ref, value = self.cache[max(0, idx - 1)]
        value += self._integrate(ref, u, self.err)
        self.cache.insert(idx, (u, value))
        return value

    def length(self, a, b):
        if b < a:
            a, b = b, a
        return self.cumulative(b) - self.cumulative(a)

    def cached_intervals(self):
        return [(a[0], b[0], b[1] - a[1]) for a, b in zip(self.cache[:-1], self.cache[1:])]


def arc_length_between(a, b, P, p, U, err=1e-8, n_start=4, n_max=256):
    return BSplineArcLength(P, p, U, err, n_start, n_max)._integrate(a, b, err)


def arc_length_bspline_gauss(P, p, U, err=1e-6, n_start=4, n_max=256):
    calc = BSplineArcLength(P, p, U, err, n_start, n_max)
    intervals = list(zip(calc.unique_knots[:-1], calc.unique_knots[1:]))
    total, used, errors = 0., [], []
    for a, b in intervals:
        n = n_start
        prev = gauss_legendre_integral(calc.speed, a, b, n)
        while True:
            curr = gauss_legendre_integral(calc.speed, a, b, 2 * n)
            diff = abs(curr - prev)
            if diff < err / len(intervals) or 2 * n >= n_max:
                total += curr
                used.append(2 * n)
                errors.append(diff)
                break
            prev, n = curr, 2 * n
    return total, dict(intervals=intervals, n_used=used, local_err=errors)


def find_u_by_arc_length(u0, L_target, P, p, U, err=1e-6, max_iter=100,
                         inner_err_ratio=.01, verbose=False):
    """Original Newton then bisection algorithm, including its query cache."""
    if L_target <= 0:
        return u0, de_boor(u0, P, p, U), 0., dict(method='trivial', iter=0, abs_err=0.)
    end = float(U[-1])
    if u0 >= end - 1e-15:
        raise ValueError("start is at the endpoint")
    calc = BSplineArcLength(P, p, U, max(err * inner_err_ratio, 1e-14))
    total = calc.length(u0, end)
    if L_target > total + err:
        raise ValueError("unreachable target")
    u = min(max(u0 + (end - u0) * L_target / total, u0 + 1e-10), end - 1e-10)
    length = calc.length(u0, u)
    method = 'newton'
    for it in range(max_iter):
        f = length - L_target
        if abs(f) <= err:
            break
        speed = calc.speed(u)
        if speed < 1e-12:
            method = 'bisection'
            break
        next_u = u - f / speed
        if not u0 < next_u < end or not np.isfinite(next_u):
            method = 'bisection'
            break
        next_length = calc.length(u0, next_u)
        if abs(next_length - L_target) >= .99 * abs(f):
            method = 'bisection'
            break
        u, length = next_u, next_length
    if method == 'newton' and abs(length - L_target) <= err:
        return u, de_boor(u, P, p, U), length, dict(method=method, iter=it+1, abs_err=abs(length-L_target))
    lo, hi = u0, end
    for it in range(max_iter):
        u = (lo + hi) / 2
        length = calc.length(u0, u)
        if abs(length - L_target) <= err:
            break
        if length < L_target:
            lo = u
        else:
            hi = u
        if hi - lo < 1e-15:
            break
    return u, de_boor(u, P, p, U), length, dict(method='bisection', iter=it+1, abs_err=abs(length-L_target))


if __name__ == '__main__':
    d, e, P, U, geo, curve = fit_bspline_by_bisection(
        (500., -150.), (0., 0.), 150., .04, .05, branch=-1, d_init=20.)
    length, info = arc_length_bspline_gauss(P, 4, U)
    u, pt, step, inv = find_u_by_arc_length(.5, 1.2, P, 4, U)
    print(f'd={d:.12g}, controls={len(P)}, max_error={e:.12g}')
    print(f'length={length:.12g}, GL difference estimate={sum(info["local_err"]):.3e}')
    print(f'u1={u:.15g}, point={pt}, length={step:.12g}, residual={inv["abs_err"]:.3e}')

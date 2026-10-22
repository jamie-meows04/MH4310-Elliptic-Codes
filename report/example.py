r"""
example_twists_sage.py -- SageMath version of example_twists.py.

The F_251 example of Construction 18(a) and Theorem 17(d),(e):
  E1 : y^2 = x^3 + x + 100 over F_251,  N = 280,  T = (65, 14) of order 5,
  phi : E1 -> E2 the Velu isogeny with kernel <T>,  H = 2 E1(F_251) (index 2),
  D = the 28 points of E2 whose fibres lie outside H,  G = k(O),  k odd.
Checks:
  * C' = C_L(phi^*D, phi^*G) is MDS (Theorem 4) and quasi-cyclic under T;
  * the twists C_L(D, G_j), G_0 = k(O), G_j = (k-1)(O) + ([j]R), are MDS
    (Theorem 4 and the direct "every k columns independent" test);
  * the zeta^j-eigenspace of translation by T on C', restricted to one
    coordinate per fibre, is a coordinate scaling of C_L(D, G_j).

Generator matrices are built from explicit Riemann-Roch bases:
  L(m O)          = < x^i (2i <= m),  x^i y (2i + 3 <= m) >,
  L((m-1)O + (S)) = L((m-1) O) + < (y + y_S)/(x - x_S) >     (S != O),
and, with psi the kernel polynomial of phi (div psi = sum_{T in K}(T) - ell(O)
for ell odd),  L(phi^*(k O)) = psi^{-k} L(ell k O).  So the rows
f(P)/psi(P)^k generate C' exactly, with no scaling.

How to run:
  * SageMath:          sage -python example_twists_sage.py [k]
  * pip "passagemath": pip install passagemath-schemes passagemath-singular
                       python example_twists_sage.py [k]
k must be odd (default 3).
"""
import sys
import itertools
from collections import defaultdict

try:
    from sage.all import *                       # full SageMath
except ImportError:                              # pip-installed passagemath
    from sage.all__sagemath_schemes import *
    from sage.all__sagemath_singular import *

p, A, B, ell = 251, 1, 100, 5
F = GF(p)


# ---------------------------------------------------------------------------
# Riemann-Roch bases as functions of a rational affine point (x, y)
# ---------------------------------------------------------------------------

def rr_functions(m, S=None):
    """Functions P -> F spanning L(m O) (S = None) or L((m-1) O + (S))."""
    if S is not None and m == 1:                 # L((S)) = constants
        return [lambda P: F(1)]
    mm = m if S is None else m - 1
    fs = [lambda P, i=i: P[0] ** i for i in range(mm // 2 + 1)]
    fs += [lambda P, i=i: P[0] ** i * P[1] for i in range((mm - 3) // 2 + 1)]
    if S is not None:
        xs, ys = S.xy()
        fs.append(lambda P: (P[1] + ys) / (P[0] - xs))
    assert len(fs) == m
    return fs


def gen_matrix(points, m, S=None, weight=None):
    """Rows (w(P) f(P))_{P in points} for f in the basis above."""
    xy = [P.xy() for P in points]
    w = [F(1)] * len(points) if weight is None else [weight(P) for P in xy]
    return matrix(F, [[w[i] * f(xy[i]) for i in range(len(xy))]
                      for f in rr_functions(m, S)])


# ---------------------------------------------------------------------------
# MDS tests
# ---------------------------------------------------------------------------

def subset_sum_count(D_pts, S, k):
    """Theorem 4: number of k-subsets of D_pts summing to S (0 <=> MDS)."""
    O = S.curve()(0)
    dp = [defaultdict(int) for _ in range(k + 1)]
    dp[0][O] = 1
    for P in D_pts:
        for j in range(k, 0, -1):
            for T, c in list(dp[j - 1].items()):
                dp[j][T + P] += c
    return dp[k][S]


def every_k_columns_independent(M, max_sets=200000):
    """Direct MDS test: every k columns of the k x n matrix M independent.
    Returns None if there are more than max_sets column sets."""
    k, n = M.nrows(), M.ncols()
    if binomial(n, k) > max_sets:
        return None
    return all(M.matrix_from_columns(list(c)).rank() == k
               for c in itertools.combinations(range(n), k))


def is_scaling(GA, GB, tries=50):
    """Is there v with all v_i != 0 and v * rowspace(GA) = rowspace(GB)?
    v * A lies in B iff sum_i h_i g_i v_i = 0 for every row g of GA and
    every parity check h of B: linear equations in v."""
    if GA.rank() != GB.rank():
        return False
    n = GA.ncols()
    HB = GB.right_kernel().basis_matrix()
    rows = [[h[i] * g[i] for i in range(n)] for g in GA.rows() for h in HB.rows()]
    Kv = matrix(F, rows).right_kernel() if rows else (F ** n)
    if Kv.dimension() == 0:
        return False
    for _ in range(tries):
        v = Kv.random_element()
        if all(c != 0 for c in v):
            return True
    return False


# ---------------------------------------------------------------------------
# The example
# ---------------------------------------------------------------------------

def main(k=3):
    assert k % 2 == 1, "Construction 18(a) needs k odd"
    E1 = EllipticCurve(F, [A, B])
    N = E1.cardinality()
    T = E1(65, 14)
    assert T.order() == ell
    # Pin the codomain to the report's model of E2, so that point
    # coordinates on E2 agree with the report.  Sage composes the Velu
    # isogeny with an isomorphism onto this curve if necessary.
    E2 = EllipticCurve(F, [231, 13])
    phi = E1.isogeny(T, codomain=E2)
    assert phi.codomain() == E2
    O1, O2 = E1(0), E2(0)
    print("E1:", E1, "  N =", N, "  group:", E1.abelian_group())
    print("E2:", E2)

    # index-2 subgroup H = 2 E1(F_p); K lies in H; fibres are K-cosets
    pts1 = E1.points()
    H = set(2 * P for P in pts1)
    assert 2 * len(H) == N, "2E1(F_p) does not have index 2"
    K = [P for P in pts1 if phi(P).is_zero()]
    assert len(K) == ell and all(P in H for P in K)
    fib = defaultdict(list)
    for P in pts1:
        fib[phi(P)].append(P)
    assert all(len(Fb) == ell and len({P in H for P in Fb}) == 1 for Fb in fib.values())

    # D: points of E2 whose fibre lies outside H; fibres ordered P~, P~+T, ...
    D2 = sorted(Q for Q, Fb in fib.items() if not Q.is_zero() and Fb[0] not in H)
    reps = [min(fib[Q]) for Q in D2]
    D1 = [r + m * T for r in reps for m in range(ell)]
    n2, n1, k1 = len(D2), len(D1), ell * k
    print("D: %d points of E2" % n2, [Q.xy() for Q in D2])

    # ker(phi-hat) and zeta
    phihat = phi.dual()                          # E2 -> E1 (up to isomorphism)
    if phihat.codomain() != E1:                  # make the codomain literally E1
        phihat = phihat.codomain().isomorphism_to(E1) * phihat
    kerhat = [Q for Q in E2.points() if not Q.is_zero() and phihat(Q).is_zero()]
    assert len(kerhat) == ell - 1, "ker(phi-hat) is not rational"
    zeta = F(20)
    assert zeta.multiplicative_order() == ell

    # C' exactly: rows f(P) / psi(x(P))^k, f in L(5k O)
    psi = phi.kernel_polynomial()
    Mc = gen_matrix(D1, k1, weight=lambda P: psi(P[0]) ** (-k))
    assert Mc.rank() == k1
    cnt = subset_sum_count(D1, O1, k1)           # Sigma(phi^*(kO)) = [k] sum(K) = O
    print("\nC' = [%d, %d, %d]  (Theorem 4 count %d => %s)"
          % (n1, k1, n1 - k1 + (cnt == 0), cnt, "MDS" if cnt == 0 else "not MDS"))
    sig = [i * ell + (m + 1) % ell for i in range(n2) for m in range(ell)]
    Mt = Mc.matrix_from_columns(sig)             # translation by T
    print("C' invariant under translation by T (quasi-cyclic):",
          Mc.stack(Mt).rank() == k1)

    # eigen-subcodes, restricted to the coordinates P~_i
    eig = []
    for e in range(ell):
        lam = zeta ** e
        X = (Mt - lam * Mc).left_kernel().basis_matrix()
        W = X * Mc
        eig.append(W.matrix_from_columns(list(range(0, n1, ell))))

    # choose the generator R of ker(phi-hat) so that zeta^1 <-> [1]R
    R0 = kerhat[0]
    twist = lambda S: gen_matrix(D2, k, None if S.is_zero() else S)
    R = -next(m * R0 for m in range(1, ell) if is_scaling(twist(m * R0), eig[1]))
    print("R =", R.xy())

    print()
    for j in range(ell):
        S = j * R
        Gj = "%d(O)" % k if S.is_zero() else "%d(O) + (%s, %s)" % ((k - 1,) + S.xy())
        Mj = twist(S)
        c = subset_sum_count(D2, S, k)
        print("twist j=%d: G_j = %-20s [%d, %d, %d]  Thm 4: %-7s direct column test: %-7s"
              " zeta^%d-eigenspace of C' is a scaling of it: %s"
              % (j, Gj, n2, k, n2 - k + (c == 0), "MDS" if c == 0 else "not MDS",
                 {True: "MDS", False: "not MDS", None: "skipped"}[every_k_columns_independent(Mj)],
                 j, eig[j].rank() == k and is_scaling(Mj, eig[j])))


if __name__ == "__main__":
    main(15)
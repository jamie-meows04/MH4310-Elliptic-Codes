"""
isogeny_codes.py -- testing the four isogeny transformations of elliptic codes
================================================================================

MH4310 project, Section 7.  For an isogeny  phi : E1 -> E2  we test

  (1)  C_L(phi^* D, phi^* G)        D, G on E2   (pull back the divisors)
  (2)  C_L(phi_* D, phi_* G)        D, G on E1   (push forward the divisors)
  (3)  phi^* C_L(D, G)              D, G on E2   (pull back the functions)
  (4)  phi_* C_L(D, G)              D, G on E1   (push forward the functions)

and, for method (1), the extra structure that comes from ker(phi) acting on E1
by translations (see method_1_structure below).

How to run
----------
  * SageMath:          sage -python isogeny_codes.py
  * Sage notebook:     paste the file into a cell (it is plain Python, no ^).
  * pip "passagemath": pip install passagemath-schemes passagemath-singular
                       python isogeny_codes.py
Tested on passagemath 10.8 (pip Sage) with Python 3.11.  Runtime about a minute.

What Sage does and does not do
------------------------------
EllipticCurve objects know points, the group law and isogenies (Velu's
formulae), but they have NO Riemann-Roch spaces.  Riemann-Roch lives in Sage's
function-field module, which works for any curve over any finite field
(Hess's algorithm).  So for each curve we build its function field
        F_q(x)[y] / (y^2 + a1 x y + a3 y - x^3 - a2 x^2 - a4 x - a6)
and keep a dictionary translating between elliptic-curve points and the
degree-1 places of that function field.  (A full Sage install also has
codes.EvaluationAGCode(places, G), which can cross-check ag_code() below.)

Minimum distances: a full Sage install computes them with GAP via
C.minimum_distance(), but that is slow for long codes.  We instead use the
subset-sum criterion of our report (Theorem 4): an elliptic code C_L(D,G) has
d = n-k if some k of the points of D add up to Sigma(G), and d = n-k+1 (MDS)
otherwise.  This is exact for every code below that is built from D and G.

Notation in comments
--------------------
  ell = deg(phi),  K = ker(phi) (a subgroup of E1 of order ell),
  "fibre over Q" = phi^{-1}(Q) = the ell points of E1 that phi sends to Q.
For separable phi with K rational, each fibre is a coset  P~ + K.
"""

import itertools
from collections import defaultdict

try:
    from sage.all import *                       # full SageMath
except ImportError:                              # pip-installed passagemath
    from sage.all__sagemath_schemes import *
    from sage.all__sagemath_singular import *
from sage.coding.linear_code import LinearCode

# `from sage.all import *` defines its own function called `random`, which
# hides Python's random module.  So we import Python's module under another
# name AFTER the Sage import.
import random as pyrandom

pyrandom.seed(4310)          # fixed seed, so every run prints the same numbers


# ---------------------------------------------------------------------------
# 1.  Elliptic curve  <->  function field
# ---------------------------------------------------------------------------

class EllipticFF:
    """The function field of an elliptic curve E in Weierstrass form, together
    with the dictionary  rational point <-> degree-1 place.

    Why we need this: Sage's EllipticCurve knows the group law but not
    Riemann-Roch spaces, while Sage's FunctionField knows Riemann-Roch spaces
    but not the group law.  This class lets us use both on the same curve."""

    def __init__(self, E, name='x'):
        self.E = E
        self.F = E.base_field()                  # the finite field F_q
        a1, a2, a3, a4, a6 = E.a_invariants()

        # Step 1: the rational function field F_q(x).
        Kx = FunctionField(self.F, name)
        X = Kx.gen()

        # Step 2: adjoin y satisfying the Weierstrass equation.  L is then the
        # function field F_q(E) = F_q(x, y) of the report's Section 2.
        Y = PolynomialRing(Kx, 'Y').gen()
        self.L = Kx.extension(Y**2 + a1*X*Y + a3*Y - (X**3 + a2*X**2 + a4*X + a6), 'y')

        # IMPORTANT: X lives in the smaller field Kx = F_q(x).  We must coerce
        # it into L before using it, otherwise Sage computes divisors on the
        # projective LINE rather than on the elliptic curve.
        self.x = self.L(X)
        self.y = self.L.gen()

        # Step 3: match points with places.  A place of degree 1 is the same
        # thing as an F_q-rational point, and evaluating x and y at it gives
        # the point's coordinates.  The point at infinity is the unique
        # infinite place.
        self.inf_place = self.L.places_infinite()[0]
        self._place = {E(0): self.inf_place}     # dictionary: point -> place
        for p in self.L.places_finite():
            if p.degree() == 1:
                P = E(self.x.evaluate(p), self.y.evaluate(p))
                self._place[P] = p
        # Sanity check: we found every rational point exactly once.
        assert len(self._place) == E.cardinality()

    def place(self, P):
        """The place of L corresponding to the rational point P of E."""
        return self._place[self.E(P)]

    def divisor(self, pts):
        """Turn a divisor written as a dictionary {point: multiplicity}, which
        is the natural way to write G = sum n_P (P), into a Sage divisor."""
        D = self.L.divisor_group().zero()
        for P, m in pts.items():
            D += self.place(P).divisor(Integer(m))   # adds m * (P)
        return D

    def riemann_roch_basis(self, G):
        """A basis of the Riemann-Roch space L(G).  G is a Sage divisor, and
        places of any degree are allowed (useful when phi^*G has
        non-rational points)."""
        return list(G.basis_function_space())


def ag_code(EF, D_pts, G):
    """The evaluation code C_L(D, G).

    D_pts : list of distinct rational points P_1, ..., P_n (in this order)
    G     : Sage divisor on EF.L whose support avoids the P_i
    Returns (the code, the basis f_1, ..., f_k of L(G) used to build it).
    Row j of the generator matrix is f_j(D) = (f_j(P_1), ..., f_j(P_n))."""
    places = [EF.place(P) for P in D_pts]
    assert len(set(places)) == len(places), "D must be a sum of DISTINCT points"
    for p in places:
        assert G.multiplicity(p) == 0, "supp D and supp G must be disjoint"
    basis = EF.riemann_roch_basis(G)
    M = matrix(EF.F, [[f.evaluate(p) for p in places] for f in basis])
    return LinearCode(M), basis


def divisor_sum(EF, G):
    """Sigma(G) = the sum of the points of G under the group law, counted with
    multiplicity.  Only works when every point of G is rational."""
    inv = {p: P for P, p in EF._place.items()}   # reverse dictionary: place -> point
    S = EF.E(0)
    for p, m in G.list():                        # G.list() = [(place, multiplicity), ...]
        assert p.degree() == 1, "Sigma() here needs a rational support"
        S += Integer(m) * inv[p]                 # m*P is [m]P in Sage
    return S


# ---------------------------------------------------------------------------
# 2.  Minimum distance
# ---------------------------------------------------------------------------

def elliptic_distance(D_pts, SigmaG, k):
    """Exact minimum distance of an elliptic code via Theorem 4.

    Returns (d, count), where count = number of k-element subsets of supp D
    whose points add up to Sigma(G).  If count > 0 then d = n-k, otherwise
    d = n-k+1 (MDS).  Also A_{n-k} = (q-1)*count.

    Why not just try every subset?  For n = 36, k = 9 there are about 94
    million subsets.  Instead we use DYNAMIC PROGRAMMING: we read the points
    one at a time and keep a table

        dp[j][S] = number of ways to pick j of the points read so far
                   so that they add up to S.

    When a new point P arrives, every old way of picking j-1 points with sum S
    extends to a way of picking j points with sum S + P.  So
        dp[j][S + P] += dp[j-1][S].
    At the end dp[k][Sigma(G)] is the answer.  The table has at most (k+1)*N
    entries (N = #E(F_q)), so the work is about n*k*N steps, a few thousand
    here instead of millions."""
    n = len(D_pts)
    O = SigmaG.curve()(0)                        # identity of the group
    # dp[j] is a dictionary {group element S: count}.  defaultdict(int)
    # returns 0 for missing keys, so we never need to check "is S in dp[j]?".
    dp = [defaultdict(int) for _ in range(k + 1)]
    dp[0][O] = 1                                 # one way to pick 0 points: sum = O
    for P in D_pts:
        # Go DOWNWARDS in j, so P is used at most once.  (Going upwards would
        # let dp[j] be updated from a dp[j-1] that already used P.)
        for j in range(k, 0, -1):
            # list(...) takes a snapshot, since we must not change a
            # dictionary while looping over it.
            for S, c in list(dp[j - 1].items()):
                dp[j][S + P] += c
    count = dp[k][SigmaG]
    return (n - k if count > 0 else n - k + 1), count


def brute_distance(C, limit=10**6):
    """Minimum distance by listing every codeword.  Only for small codes
    (q^k <= limit); beyond that we try Sage/GAP and give up (None) if GAP is
    not installed.  We use it to cross-check elliptic_distance."""
    q, k = C.base_field().cardinality(), C.dimension()
    if q**k > limit:
        try:
            return C.minimum_distance()
        except Exception:
            return None
    return min(v.hamming_weight() for v in C if not v.is_zero())


def params(C, d=None):
    """Pretty-print [n, k, d]."""
    return "[%d, %d, %s]" % (C.length(), C.dimension(), "?" if d is None else d)


# ---------------------------------------------------------------------------
# 3.  Equivalence by coordinate scaling (detects iso-duality etc.)
# ---------------------------------------------------------------------------

def scaling_to(A, B, tries=60):
    """Look for v with every v_i != 0 and  v * A = B, where
    v * A = {(v_1 a_1, ..., v_n a_n) : a in A}.  Returns v or None.

    Idea: v * A is contained in B  iff  every row g of A's generator matrix,
    scaled by v, is orthogonal to every row h of B's parity-check matrix, i.e.
        sum_i h_i g_i v_i = 0.
    These are LINEAR equations in the unknowns v_1, ..., v_n, so the solutions
    form a vector space (a matrix kernel).  We still need a solution with no
    zero entry; a random element of the kernel almost always works, so we try
    a few.  Containment plus equal dimension gives v * A = B."""
    if A.length() != B.length() or A.dimension() != B.dimension():
        return None
    F, n = A.base_field(), A.length()
    GA, HB = A.generator_matrix(), B.parity_check_matrix()
    # One equation (one row) for each pair (g, h); column i holds h_i * g_i.
    rows = [[h[i] * g[i] for i in range(n)] for g in GA.rows() for h in HB.rows()]
    K = matrix(F, rows).right_kernel() if rows else (F**n)
    if K.dimension() == 0:
        return None
    basis = K.basis()
    for _ in range(tries):
        v = sum(F.random_element() * b for b in basis)   # random kernel element
        if all(c != 0 for c in v):
            return v
    return None


# ---------------------------------------------------------------------------
# 4.  Pullbacks and pushforwards along an isogeny phi : E1 -> E2
# ---------------------------------------------------------------------------

class IsogenyTools:
    """Everything we need to move functions, points and divisors along phi."""

    def __init__(self, phi):
        self.phi = phi
        self.E1, self.E2 = phi.domain(), phi.codomain()
        self.EF1, self.EF2 = EllipticFF(self.E1, 'x1'), EllipticFF(self.E2, 'x2')
        self.deg = phi.degree()                       # ell

        # Sage describes phi by two rational functions:
        #     phi(x, y) = (X(x, y), Y(x, y)).
        # They come as fractions of polynomials in x, y.  We substitute the
        # function-field elements x, y of E1 into numerator and denominator,
        # which turns X and Y into elements of F_q(E1).  These are exactly the
        # pullbacks phi^*(x2) and phi^*(y2) of the coordinate functions of E2.
        Xr, Yr = phi.rational_maps()
        sub = lambda r: r.numerator()(self.EF1.x, self.EF1.y) / r.denominator()(self.EF1.x, self.EF1.y)
        self.X, self.Y = sub(Xr), sub(Yr)

        # Fibres: for every rational point Q of E2, the list of rational
        # points P of E1 with phi(P) = Q.  We build this dictionary in ONE pass
        # over E1(F_q) (fine for small q).  If ker(phi) is rational, every
        # fibre in the image has exactly ell points and the others are empty.
        self.pre = defaultdict(list)
        for P in self.E1.points():
            self.pre[phi(P)].append(P)
        self.image = sorted(self.pre.keys())          # phi(E1(F_q)), a subgroup
        self.kernel = self.pre[self.E2(0)]            # the fibre over infinity

    # -- functions -----------------------------------------------------------
    def pull_function(self, f):
        """phi^* f = f o phi = f(X, Y), for f in F_q(E2).

        Sage stores f as c_0(x) + c_1(x) y with c_0, c_1 in F_q(x)
        (f.list() returns [c_0, c_1]).  We substitute x -> X, y -> Y."""
        ev = lambda c: c.numerator()(self.X) / c.denominator()(self.X)
        return sum(ev(c) * self.Y**i for i, c in enumerate(f.list()))

    def push_value(self, f, Q):
        """(phi_* f)(Q) = product of f over the fibre phi^{-1}(Q).

        This is the norm of f evaluated at Q (separable case).  Needs the whole
        fibre to be rational, and f to have no pole on it."""
        P_list = self.pre[Q]
        assert len(P_list) == self.deg, "preimages of Q are not all rational"
        return prod(f.evaluate(self.EF1.place(P)) for P in P_list)

    # -- divisors ------------------------------------------------------------
    def pull_point(self, Q):
        """phi^*(Q) as a Sage divisor on E1.

        This works for ANY isogeny kernel: when the fibre is not rational,
        Sage simply returns places of degree > 1.  It works by finding zeros
        and poles, never by listing points:

        * Q = infinity: x2 has a double pole at infinity and no other pole, so
          phi^*(x2) = X has a double pole exactly at the fibre over infinity.
          Halving the pole divisor of X gives phi^*(infinity).
        * Q = (a, b): x2 - a vanishes at Q and -Q, and y2 - b vanishes at Q
          and two other points.  Their only COMMON zero is Q (with
          multiplicity 1), so the common zeros of X - a and Y - b, taking the
          smaller multiplicity at each place, form phi^*(Q).  (This uses that
          phi is separable, i.e. unramified.)"""
        Q = self.E2(Q)
        if Q.is_zero():
            Pd = self.X.divisor_of_poles()
            return sum((p.divisor(Integer(m // 2)) for p, m in Pd.list()),
                       self.EF1.L.divisor_group().zero())
        a, b = Q.xy()
        Z1 = (self.X - a).divisor_of_zeros()
        Z2 = (self.Y - b).divisor_of_zeros()
        out = self.EF1.L.divisor_group().zero()
        for p in Z1.support():
            m = min(Z1.multiplicity(p), Z2.multiplicity(p))
            if m > 0:
                out += p.divisor(m)
        return out

    def pull_divisor(self, pts):
        """phi^* of a divisor {point: multiplicity} on E2.  Pullback is
        additive, so phi^*(sum n_Q (Q)) = sum n_Q phi^*(Q)."""
        return sum((Integer(m) * self.pull_point(Q) for Q, m in pts.items()),
                   self.EF1.L.divisor_group().zero())

    def pull_points(self, D_pts):
        """phi^* D as a LIST of rational points (needed for evaluation).

        The list is ordered fibre by fibre: the ell points over D_pts[0],
        then the ell points over D_pts[1], and so on.  Raises an error if some
        fibre is not rational, since then phi^*D is not a valid evaluation
        divisor."""
        out = []
        for Q in D_pts:
            P_list = self.pre[Q]
            if len(P_list) != self.deg:
                raise ValueError("phi^*(%s) is not a sum of %d rational points" % (Q, self.deg))
            out += P_list
        return out

    def push_points(self, pts):
        """phi_* of a divisor {point: multiplicity} on E1: just apply phi to
        each point and add up multiplicities that land on the same point."""
        out = defaultdict(int)
        for P, m in pts.items():
            out[self.phi(P)] += m
        return dict(out)


# ---------------------------------------------------------------------------
# 5.  The experiments
# ---------------------------------------------------------------------------

def banner(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


def setup(q=31, a=(0, 0, 0, 0, 5), ell=3):
    """Pick E1 : y^2 = x^3 + 5 over F_31 (39 points, cyclic) and a rational
    point T of order 3.  Velu's formulae (E1.isogeny(T)) give the isogeny
    phi : E1 -> E2 with kernel {O, T, 2T}.  Because T is rational, every
    rational point of phi(E1(F_q)) has a fully rational fibre."""
    F = GF(q)
    E1 = EllipticCurve(F, list(a))
    T = next(P for P in E1.points() if P.order() == ell)   # first point of order ell
    phi = E1.isogeny(T)
    I = IsogenyTools(phi)
    banner("SET-UP")
    print("E1 :", E1, "   #E1 =", E1.cardinality(), E1.abelian_group())
    # Isogenous curves over F_q have the same number of points (Tate), so
    # #E2 = #E1 always.  phi(E1(F_q)) is a subgroup of index ell in E2(F_q).
    print("E2 :", I.E2, "   #E2 =", I.E2.cardinality())
    print("phi:", phi, "\n     x ->", phi.rational_maps()[0])
    print("ker phi =", I.kernel, " (rational, so preimages of image points are rational)")
    print("#phi(E1(F_q)) =", len(I.image), "  (index", I.E2.cardinality() // len(I.image), "in E2(F_q))")
    return I


def method_1_and_3(I, k=3, Gpts=None, D2=None, show=True):
    """(1) C_L(phi^*D, phi^*G) and (3) phi^*C_L(D,G), for D, G on E2.

    Defaults: G = k(infinity) and D = every affine point of phi(E1(F_q)),
    the largest D whose pullback is rational."""
    E2, EF2, EF1 = I.E2, I.EF2, I.EF1
    O2 = E2(0)
    Gpts = Gpts or {O2: k}
    G2 = EF2.divisor(Gpts)
    if D2 is None:
        D2 = [Q for Q in I.image if Q not in Gpts]
    # The original code C on E2 and its exact distance.
    C, basis2 = ag_code(EF2, D2, G2)
    d, cnt = elliptic_distance(D2, divisor_sum(EF2, G2), k)

    # (1) Pull back the divisors.  phi^*D has ell*n points and phi^*G has
    # degree ell*k, so C1 is an [ell n, ell k, >= ell(n - k)] code.
    D1 = I.pull_points(D2)
    G1 = I.pull_divisor(Gpts)
    C1, _ = ag_code(EF1, D1, G1)
    # Theorem 4 needs Sigma(phi^*G), which needs rational points in phi^*G.
    rationalG1 = all(p.degree() == 1 for p in G1.support())
    d1 = elliptic_distance(D1, divisor_sum(EF1, G1), C1.dimension())[0] if rationalG1 else brute_distance(C1)

    # (3) Pull back the functions: row j is (phi^* f_j)(phi^* D) with f_j in
    # L(G).  Since (phi^* f)(P) = f(phi(P)), every fibre just repeats the
    # value f(Q), so we expect the ell-fold repetition of C.
    M3 = matrix(EF1.F, [[I.pull_function(f).evaluate(EF1.place(P)) for P in D1] for f in basis2])
    C3 = LinearCode(M3)
    d3 = brute_distance(C3)

    if show:
        banner("METHODS (1) AND (3): D, G on E2, pulled back to E1")
        print("original  C  = C_L(D,G) on E2      :", params(C, d), " MDS" if d == C.length() - k + 1 else " near-MDS",
              "  (brute force d = %s)" % brute_distance(C))
        print("(1)       C1 = C_L(phi^*D, phi^*G) :", params(C1, d1),
              " MDS" if d1 == C1.length() - C1.dimension() + 1 else " near-MDS")
        print("(3)       C3 = phi^*C_L(D,G)       :", params(C3, d3))
        # Check the repetition claim column by column.  D1 is ordered fibre by
        # fibre, so columns ell*i, ..., ell*i + ell - 1 of M3 should all equal
        # column i of C's generator matrix M.
        M = matrix(EF2.F, [[f.evaluate(EF2.place(Q)) for Q in D2] for f in basis2])
        rep = all(M3.column(I.deg * i + t) == M.column(i)
                  for i in range(len(D2)) for t in range(I.deg))
        print("    C3 is the %d-fold repetition of C       :" % I.deg, rep)
        # phi^* L(G) is inside L(phi^* G), so C3 should sit inside C1.
        print("    C3 is a subcode of C1                    :", all(r in C1 for r in C3.generator_matrix().rows()))

        # Duality commutes with pullback.  On an elliptic curve
        # C_L(D,G)^perp is a scaling of C_L(D,H') for ANY H' of degree n-k with
        # Sigma(H') = Sigma(D) - Sigma(G).  We take the two-point normal form
        # H' = (n-k-1)(infinity) + (Sigma(D) - Sigma(G)) and check that
        # C1^perp is likewise a scaling of C_L(phi^*D, phi^*H').
        SigmaH = sum(D2, O2) - divisor_sum(EF2, G2)
        Hpts = {O2: len(D2) - k - 1}
        Hpts[SigmaH] = Hpts.get(SigmaH, 0) + 1
        if all(Q not in D2 for Q in Hpts):
            Hcode, _ = ag_code(EF2, D2, EF2.divisor(Hpts))
            H1code, _ = ag_code(EF1, D1, I.pull_divisor(Hpts))
            print("    C^perp  ~ C_L(D, H')        by scaling  :", scaling_to(Hcode, C.dual_code()) is not None)
            print("    C1^perp ~ C_L(phi^*D, phi^*H') by scaling:", scaling_to(H1code, C1.dual_code()) is not None)
        # Iso-duality only makes sense when n = 2k.
        if C.length() == 2 * C.dimension():
            print("    C  iso-dual (by scaling)                 :", scaling_to(C, C.dual_code()) is not None)
            print("    C1 iso-dual (by scaling)                 :", scaling_to(C1, C1.dual_code()) is not None)
            # Hull dimension = k - rank(G G^T) for a generator matrix G.
            hull = lambda Z: Z.dimension() - (Z.generator_matrix() * Z.generator_matrix().transpose()).rank()
            print("    hull dimensions  dim Hull(C), dim Hull(C1):", hull(C), hull(C1))
    return C, d, C1, d1


def method_1_trials(I, trials=150):
    """Random D, G on E2.  Prediction: if C1 is MDS then C is MDS.

    Reason: if k points S of D add up to Sigma(G), the ell*k points over S add
    up to Sigma(phi^*G), so C1 fails Theorem 4 whenever C does."""
    banner("METHOD (1): %d random trials  (prediction: C1 MDS  =>  C MDS)" % trials)
    tally = defaultdict(int)          # counts each (C MDS?, C1 MDS?) outcome
    img = list(I.image)
    for _ in range(trials):
        k = pyrandom.randint(2, 5)
        # G = (k-1)(Q0) + (Q1), the two-point normal form, with Q0, Q1 in the
        # image so that phi^*G is rational and Theorem 4 applies upstairs.
        Q0, Q1 = pyrandom.choice(img), pyrandom.choice(img)
        Gpts = {Q0: k - 1}
        Gpts[Q1] = Gpts.get(Q1, 0) + 1
        pool = [Q for Q in img if Q not in Gpts]
        n = pyrandom.randint(k + 1, len(pool))
        D2 = pyrandom.sample(pool, n)          # n distinct points, no repeats
        SG = (k - 1) * Q0 + Q1                 # Sigma(G)
        d = elliptic_distance(D2, SG, k)[0]
        D1 = I.pull_points(D2)
        # Sigma(phi^*G): add up every point of every fibre, with multiplicity.
        SG1 = sum((Integer(m) * P for Q, m in Gpts.items() for P in I.pre[Q]), I.E1(0))
        d1 = elliptic_distance(D1, SG1, I.deg * k)[0]
        tally[(d == n - k + 1, d1 == I.deg * (n - k) + 1)] += 1
    for (a, b), c in sorted(tally.items()):
        print("  C MDS = %-5s  C1 MDS = %-5s : %d" % (a, b, c))
    print("  counterexamples to the prediction:", tally[(False, True)])


def method_1_structure(I, k=3):
    """The extra structure of C1 = C_L(phi^*D, phi^*G), coming from K = ker phi.

    (a) For T in K we have phi(P + T) = phi(P), so translation by T maps
        phi^*D and phi^*G to THEMSELVES (not just to equivalent divisors).
        Hence T acts on C1 by a pure coordinate permutation, with no scaling,
        moving each coordinate around its fibre.  With K cyclic of order ell,
        C1 is a quasi-cyclic code of index n and co-index ell.
    (b) The codewords fixed by every T are exactly the evaluations of
        K-invariant functions, i.e. of functions pulled back from E2.  So the
        invariant subcode should be C3 = phi^*C, the repetition of C.
    (c) If ell divides q - 1, the permutation splits C1 into ell eigen-subcodes
        (one for each ell-th root of unity w^j).  Our prediction is that each
        one is, after keeping one coordinate per fibre, an elliptic code on E2
        with the SAME D and a "twisted" G_j, where
            Sigma(G_j) = Sigma(G) + [j]R,   R a generator of ker(phi-hat),
        and phi-hat : E2 -> E1 is the dual isogeny.  So C1 glues together the
        ell elliptic codes C_L(D, G_j), j = 0, ..., ell-1."""
    banner("METHOD (1) STRUCTURE: quasi-cyclic decomposition by ker(phi)")
    F, E1, E2, EF1, EF2 = I.EF1.F, I.E1, I.E2, I.EF1, I.EF2
    O2 = E2(0)
    Gpts = {O2: k}
    D2 = [Q for Q in I.image if Q not in Gpts]
    D1 = I.pull_points(D2)
    C1, _ = ag_code(EF1, D1, I.pull_divisor(Gpts))
    n1 = len(D1)

    # (a) The coordinate permutation induced by a generator T of K:
    # coordinate j (the point D1[j]) is sent to the coordinate of D1[j] + T.
    T = next(P for P in I.kernel if not P.is_zero())
    index = {P: j for j, P in enumerate(D1)}            # point -> coordinate
    sigma = [index[D1[j] + T] for j in range(n1)]
    # Permutation matrix Pm with (c * Pm)_j = c_{sigma(j)}.
    Pm = matrix(F, n1, n1)
    for j in range(n1):
        Pm[sigma[j], j] = 1
    G1 = C1.generator_matrix()
    invariant = all((r * Pm) in C1 for r in G1.rows())
    print("(a) translation by T permutes C1 (no scaling needed):", invariant)

    # (b) Fixed vectors of the permutation = left kernel of (Pm - Id).
    V = G1.row_space()
    fixed = V.intersection((Pm - identity_matrix(F, n1)).left_kernel())
    print("(b) dim of the K-invariant subcode = %d   (k = %d)" % (fixed.dimension(), k))

    # (c) Eigen-subcodes, if F_q contains the ell-th roots of unity.
    ell = I.deg
    if (F.cardinality() - 1) % ell != 0:
        print("(c) skipped: ell does not divide q - 1")
        return
    w = F.zeta(ell)                                     # primitive ell-th root of 1
    # A generator R of ker(phi-hat) on E2 (rational in this example).
    phihat = I.phi.dual()
    R = next((P for P in E2.points() if not P.is_zero() and phihat(P).is_zero()), None)
    print("    ker(phi-hat) generated by R =", R)
    # The twisted codes C_L(D, G_j) on E2, with G_j in two-point normal form
    # G_j = (k-1)(infinity) + (Sigma(G) + [j]R).
    twisted = []
    for j in range(ell):
        S = divisor_sum(EF2, EF2.divisor(Gpts)) + j * R
        Gj = {O2: k} if S.is_zero() else {O2: k - 1, S: 1}
        assert all(Q not in D2 for Q in Gj), "twist point lands in D; choose another D"
        twisted.append(ag_code(EF2, D2, EF2.divisor(Gj))[0])
    # The first coordinate of each fibre: positions 0, ell, 2*ell, ...
    reps = list(range(0, n1, ell))
    for e in range(ell):
        lam = w**e
        Ve = V.intersection((Pm - lam * identity_matrix(F, n1)).left_kernel())
        # Keep only one coordinate per fibre.  This loses nothing, because an
        # eigenvector is determined by one value per fibre (the others are that
        # value times powers of lam).
        proj = matrix(F, [[v[i] for i in reps] for v in Ve.basis()])
        Ce = LinearCode(proj) if proj.nrows() > 0 else None
        matches = [j for j in range(ell)
                   if Ce is not None and scaling_to(twisted[j], Ce) is not None]
        print("(c) eigenvalue w^%d: eigen-subcode dim %d; matches C_L(D, G_j) for j in %s"
              % (e, Ve.dimension(), matches))


def method_2_and_4(I, k=3):
    """(2) C_L(phi_*D, phi_*G) and (4) phi_* C_L(D,G), for D, G on E1."""
    E1, EF1, EF2 = I.E1, I.EF1, I.EF2
    O1 = E1(0)
    banner("METHODS (2) AND (4): D, G on E1, pushed forward to E2")
    # D1: one point from each non-kernel fibre, so that phi_*D has no repeated
    # points.  (pyrandom.randrange(ell) picks which of the ell points to use.)
    D1 = [I.pre[Q][pyrandom.randrange(I.deg)] for Q in I.image if not Q.is_zero()]
    Gpts = {O1: k}
    G1 = EF1.divisor(Gpts)
    C, basis1 = ag_code(EF1, D1, G1)
    d = elliptic_distance(D1, O1, k)[0]
    # (2) phi_* D = sum of the images, phi_* G = k(infinity) again.
    D2 = [I.phi(P) for P in D1]
    G2pts = I.push_points(Gpts)
    C2, _ = ag_code(EF2, D2, EF2.divisor(G2pts))
    d2 = elliptic_distance(D2, I.E2(0), k)[0]
    print("original C  = C_L(D,G) on E1       :", params(C, d))
    print("(2)      C2 = C_L(phi_*D, phi_*G)  :", params(C2, d2))
    # What goes wrong if D contains two points of the same fibre:
    bad = I.pre[I.image[1]][:2]
    print("    phi_* of two points in one fibre =", I.push_points({P: 1 for P in bad}),
          " -> not a valid evaluation divisor")

    # (4) Push forward the functions:  the "codeword" of f is
    #       ( (phi_* f)(phi(P)) )_{P in D} = ( product of f over the fibre of P ).
    # We list ALL f in L(G) (q^k of them) and collect the resulting vectors.
    # To make this fast, first evaluate each basis function b_j once at every
    # point of every fibre; then f = sum_j c_j b_j is evaluated by arithmetic
    # on these stored numbers, with no further calls to Sage's evaluate().
    fibres = [I.pre[I.phi(P)] for P in D1]
    vals = [[[b.evaluate(EF1.place(P)) for P in fib] for fib in fibres] for b in basis1]
    #   vals[j][i][t] = b_j at the t-th point of the i-th fibre
    F = EF1.F
    S = set()                         # a set, so duplicate vectors are stored once
    # itertools.product(list(F), repeat=k) runs over every coefficient vector
    # (c_0, ..., c_{k-1}) in F^k, i.e. over every f in L(G).
    for coeffs in itertools.product(list(F), repeat=k):
        if all(c == 0 for c in coeffs):
            continue
        v = tuple(prod(sum(coeffs[j] * vals[j][i][t] for j in range(k)) for t in range(I.deg))
                  for i in range(len(D1)))
        S.add(v)
    S = [vector(F, v) for v in S]
    Sset = set(tuple(v) for v in S)
    # A linear code is closed under addition; test 300 random pairs.
    closed = all(tuple(a + b) in Sset or (a + b).is_zero()
                 for a, b in (pyrandom.sample(S, 2) for _ in range(300)))
    span = matrix(F, S).row_space()
    print("(4)      phi_*C_L(D,G): %d distinct nonzero vectors (a linear [%d,%d] code would have %d)"
          % (len(S), len(D1), k, F.cardinality()**k - 1))
    print("    closed under addition (300 random pairs) :", closed)
    print("    dimension of its linear span             :", span.dimension())
    print("    span == C2 = C_L(phi_*D, phi_*G)         :", span == C2.generator_matrix().row_space())
    print("    min weight of a nonzero vector in the set:", min(v.hamming_weight() for v in S))


def method_2_trials(I, trials=150):
    """Random D, G on E1.  Prediction: if C2 is MDS then C is MDS.

    Reason: C2 fails Theorem 4 iff some k points of D add up to an element of
    Sigma(G) + K, which is ell targets instead of the single target Sigma(G)."""
    banner("METHOD (2): %d random trials  (prediction: C2 MDS  =>  C MDS)" % trials)
    E1 = I.E1
    tally = defaultdict(int)
    reps = list(I.image)
    for _ in range(trials):
        k = pyrandom.randint(2, 5)
        P0, P1 = pyrandom.choice(E1.points()), pyrandom.choice(E1.points())
        Gpts = {P0: k - 1}
        Gpts[P1] = Gpts.get(P1, 0) + 1
        # Avoid the fibres of supp G, so that phi(supp D) misses phi(supp G).
        bad_fibres = {I.phi(P) for P in Gpts}
        pool = [Q for Q in reps if Q not in bad_fibres]
        n = pyrandom.randint(k + 1, len(pool))
        # One random point from each of n random fibres, so phi_*D is reduced.
        D1 = [pyrandom.choice(I.pre[Q]) for Q in pyrandom.sample(pool, n)]
        d = elliptic_distance(D1, (k - 1) * P0 + P1, k)[0]
        D2 = [I.phi(P) for P in D1]
        d2 = elliptic_distance(D2, I.phi((k - 1) * P0 + P1), k)[0]   # Sigma(phi_*G) = phi(Sigma G)
        tally[(d == n - k + 1, d2 == n - k + 1)] += 1
    for (a, b), c in sorted(tally.items()):
        print("  C MDS = %-5s  C2 MDS = %-5s : %d" % (a, b, c))
    print("  counterexamples to the prediction:", tally[(False, True)])


def iso_dual_check(I):
    """Prediction: the pullback of an iso-dual code is iso-dual.  We need n = 2k.
    D = the 12 affine points of phi(E1(F_q)), which add up to infinity, and
    G = 6(infinity), so Sigma(H') = Sigma(D) - Sigma(G) = infinity = Sigma(G)."""
    banner("ISO-DUALITY UNDER PULLBACK  (n = 2k on E2, deg G = k)")
    E2 = I.E2
    D2 = [Q for Q in I.image if not Q.is_zero()]
    k = len(D2) // 2
    method_1_and_3(I, k=k, Gpts={E2(0): k}, D2=D2)


def inseparable_demo(q=31):
    """Frobenius pi(x, y) = (x^p, y^p) is inseparable.  pi^*(P) is p(P), a
    single point with multiplicity p, so pi^*D is never a sum of distinct
    points and method (1) cannot even be set up."""
    banner("INSEPARABLE ISOGENY (Frobenius) : pullback is never reduced")
    E = EllipticCurve(GF(q), [0, 5])
    EF = EllipticFF(E)
    P = next(P for P in E.points() if not P.is_zero())
    a, b = P.xy()
    Z1 = (EF.x**q - a).divisor_of_zeros()
    Z2 = (EF.y**q - b).divisor_of_zeros()
    m = min(Z1.multiplicity(EF.place(P)), Z2.multiplicity(EF.place(P)))
    print("pi^*(P) has multiplicity", m, "at P  =>  pi^*D is not a sum of distinct points")


if __name__ == "__main__":
    I = setup()
    method_1_and_3(I, k=3)
    iso_dual_check(I)
    method_1_trials(I)
    method_1_structure(I, k=3)
    method_2_and_4(I, k=3)
    method_2_trials(I)
    inseparable_demo()

"""
Pullback codes over 𝔽₆₄: Example 28 of
"Isogeny Pullbacks of Elliptic Codes" (beyond length q + 1).

Coordinates of φ*D are ordered block by block: block m consists of
the points P̃ᵢ ⊕ [m]T for i = 1, …, n.
"""

from itertools import combinations

from sage.all import GF, PolynomialRing, EllipticCurve, matrix, vector

# --------------------------------------------------------------------------
# Printing helpers
# --------------------------------------------------------------------------
_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")
_SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def sup(e):
    return str(e).translate(_SUP)


def sub(e):
    return str(e).translate(_SUB)


def fe(c):
    """A field element as a power of α."""
    if c == 0:
        return "0"
    e = alog(c)
    if e == 0:
        return "1"
    return "α" if e == 1 else "α" + sup(e)


def fpt(P):
    """A point as (x, y) with coordinates as powers of α, or ∞."""
    if P.is_zero():
        return "∞"
    return "(%s, %s)" % (fe(P[0]), fe(P[1]))


def fterm(c, mono):
    """The term c·mono, e.g. α²x, xy or 1."""
    if c == 0:
        return None
    if mono == "":
        return fe(c)
    return mono if c == 1 else fe(c) + mono


def fcurve(E, name):
    """The Weierstrass equation of E."""
    a1, a2, a3, a4, a6 = E.a_invariants()
    lhs = [t for t in ("y²", fterm(a1, "xy"), fterm(a3, "y")) if t]
    rhs = [t for t in ("x³", fterm(a2, "x²"), fterm(a4, "x"), fterm(a6, ""))
           if t]
    return "%s: %s = %s" % (name, " + ".join(lhs), " + ".join(rhs))


def fpoly(f, var="α"):
    """A polynomial over 𝔽₂ in the variable var."""
    terms = []
    for e in sorted(f.exponents(), reverse=True):
        terms.append("1" if e == 0 else var + (sup(e) if e > 1 else ""))
    return " + ".join(terms)


def fvec(v, per_line=11, indent="    "):
    """A vector of field elements, several entries per line."""
    items = ["%3s" % fe(c) for c in v]
    lines = [", ".join(items[i:i + per_line])
             for i in range(0, len(items), per_line)]
    return indent + ("\n" + indent).join(lines)


def fpts(pts, per_line=4, indent="    "):
    """A list of points, several per line."""
    items = [fpt(P) for P in pts]
    lines = [", ".join(items[i:i + per_line])
             for i in range(0, len(items), per_line)]
    return indent + ("\n" + indent).join(lines)


def header(title):
    print()
    print("=" * 76)
    print(title)
    print("=" * 76)


# --------------------------------------------------------------------------
# Field and curves
# --------------------------------------------------------------------------
Fx = PolynomialRing(GF(2), 'x')
X = Fx.gen()
F = GF(2**6, 'α', modulus=X**6 + X + 1)
α = F.gen()
q = F.order()


def alog(c):
    """Exponent e with c = αᵉ."""
    return c.log(α)


def point_key(P):
    return (alog(P[0]), alog(P[1]) if P[1] != 0 else -1)


header("Example 28")
print("𝔽₆₄ = 𝔽₂(α), minimal polynomial of α:  %s = 0" % fpoly(F.modulus()))
print("    i.e. α⁶ = %s" % fpoly(F.modulus() - X**6))
# Check: α is a primitive element of 𝔽₆₄ (order 63)
assert α.multiplicative_order() == q - 1

β = α**3 + α**2 + α + 1
print("β = α³ + α² + α + 1 = %s" % fe(β))
E1 = EllipticCurve(F, [1, 0, 0, 0, β])          # E₁: y² + xy = x³ + β
pts1 = E1.points()                              # includes ∞
O1 = E1(0)
N = len(pts1)
# Check: N = #E₁(𝔽₆₄) = 72
assert N == 72
print(fcurve(E1, "E₁"))
print("N = #E₁(𝔽₆₄) = %d" % N)

E1_3 = [P for P in pts1 if 3 * P == O1]
# Check: E₁[3] ⊆ E₁(𝔽₆₄): all 9 points of order dividing 3 are rational
assert len(E1_3) == 9
print("E₁[3] ⊆ E₁(𝔽₆₄) (%d points)" % len(E1_3))

two_torsion = [P for P in pts1 if P != O1 and 2 * P == O1]
# Check: E₁ has exactly one rational point of order 2
assert len(two_torsion) == 1
print("Rational point of order 2 on E₁: %s" % fpt(two_torsion[0]))

# Kernel K = ⟨T⟩ of order ℓ = 3, and Vélu's isogeny φ: E₁ → E₂
T = E1(α**2, α)
# Check: T has order ℓ = 3
assert T.order() == 3
ℓ = 3
K = [O1, T, 2 * T]
E2 = EllipticCurve(F, [1, 0, 0, α**2, α**35])  # E₂: y² + xy = x³ + α²x + α³⁵
O2 = E2(0)
# Sage may realise φ only up to an automorphism of E₂ (here ±1). Fix the
# sign so that φ is Vélu's normalised isogeny: φ(α, α⁵⁰) = (α⁵⁴, α³⁰).
φ_sage = E1.isogeny(T, codomain=E2)
P_ref, φP_ref = E1(α, α**50), E2(α**54, α**30)
sign = 1 if φ_sage(P_ref) == φP_ref else -1


def φ(P):
    return sign * φ_sage(P)


# Check: φ is Vélu's normalised isogeny (not its negative)
assert φ(P_ref) == φP_ref
print("T = %s, of order ℓ = %d" % (fpt(T), T.order()))
print("K = ⟨T⟩ = {%s}" % ", ".join(fpt(P) for P in K))
print("φ: E₁ → E₂ is Vélu's isogeny with kernel K, where")
print("    " + fcurve(E2, "E₂"))

# --------------------------------------------------------------------------
# Eigenfunction u, divisor A = (R) − (∞), and ker φ̂
# (ker φ̂ = ⟨R⟩ decides the cases of Theorem 25(c).)
# --------------------------------------------------------------------------
ζ = α**21
# Check: ζ is a primitive ℓ-th root of unity
assert ζ.multiplicative_order() == ℓ
print("ζ = %s, of order %d" % (fe(ζ), ζ.multiplicative_order()))

xT = T[0]
u_c1, u_c0 = α**56, α**14                       # u = (y + c₁x + c₀)/(x + x(T))


def u(P):
    """u = (y + α⁵⁶x + α¹⁴)/(x + α²), for P ∉ K."""
    x, y = P[0], P[1]
    return (y + u_c1 * x + u_c0) / (x + xT)


affine_off_K = [P for P in pts1 if P not in K]
# Check: u is a ζ-eigenfunction of τ = τ_T*: u(P ⊕ T) = ζ·u(P)
assert all(u(P + T) == ζ * u(P) for P in affine_off_K if P + T not in K)
print("u = (y + %sx + %s)/(x + %s), with τu = ζu for τ = τ_T*"
      % (fe(u_c1), fe(u_c0), fe(xT)))

R = E2(α**12, α**10)
# Check: R has order ℓ
assert R.order() == ℓ
# Check: the zeros of u away from K are exactly the points over R, so div u = φ*A
assert all((u(P) == 0) == (φ(P) == R) for P in affine_off_K)
fibre_R = [P for P in pts1 if φ(P) == R]
# Check: the fibre over R consists of ℓ rational points
assert len(fibre_R) == ℓ
R2 = 2 * R
print("R = %s, and div u = φ*A with A = (R) − (∞)" % fpt(R))
print("    zeros of u: φ⁻¹(R) = {%s}" % ", ".join(fpt(P) for P in fibre_R))

ker_φ_dual = set(φ(P) for P in E1_3)
# Check: φ(E₁[3]) = ⟨R⟩
assert ker_φ_dual == {O2, R, R2}
# (Sage cannot compute φ̂ in characteristic 2, so check ⟨R⟩ = ker φ̂ directly:
#  ψ: E₂ → E₂/⟨R⟩ satisfies ψ∘φ ≅ [3], i.e. ψ = φ̂ up to isomorphism.)
ψ = E2.isogeny(R)
# Check: the isogeny ψ with kernel ⟨R⟩ maps back to (a curve isomorphic to) E₁
assert ψ.codomain().is_isomorphic(E1)
# Check: ker(ψ∘φ) ∩ E₁(𝔽₆₄) = E₁[3], so ψ∘φ ≅ [3] and ker φ̂ = ⟨R⟩
assert [P for P in pts1 if ψ(φ(P)) == ψ.codomain()(0)] == E1_3
print("[2]R = %s" % fpt(R2))
print("ker φ̂ = φ(E₁[3]) = ⟨R⟩ = {%s}"
      % ", ".join(fpt(P) for P in (O2, R, R2)))

# --------------------------------------------------------------------------
# D: the affine points of φ(E₁(𝔽_q)) other than the image of 2-torsion
# --------------------------------------------------------------------------
fibres = {}
for P in pts1:
    fibres.setdefault(φ(P), []).append(P)
image = sorted((P for P in fibres if not P.is_zero()), key=point_key)
# Check: φ(E₁(𝔽₆₄)) has N/ℓ = 24 points (23 affine points and ∞)
assert len(image) + 1 == N // ℓ
φ_two = φ(two_torsion[0])
D = [P for P in image if P != φ_two]
n = len(D)
# Check: D has n = 22 points
assert n == 22
# Check: D = −D
assert set(-P for P in D) == set(D)
print("φ(E₁(𝔽₆₄)) has %d points; φ of the 2-torsion point is %s"
      % (len(image) + 1, fpt(φ_two)))
print("D has n = %d points P₁, …, P%s (D = −D):" % (n, sub(n)))
print(fpts(D))
print("Note: R, [2]R ∈ supp D")

Pt = [min(fibres[P], key=point_key) for P in D]
print("Representatives P̃ᵢ ∈ φ⁻¹(Pᵢ):")
print(fpts(Pt))

# --------------------------------------------------------------------------
# Riemann–Roch bases
# --------------------------------------------------------------------------


def monomials(m):
    """Basis xⁱyʲ (j ≤ 1, 2i + 3j ≤ m) of 𝓛(m∞)."""
    if m < 0:
        return []
    return [(i, j) for j in (0, 1) for i in range(m // 2 + 1)
            if 2 * i + 3 * j <= m]


def ev_mono(mono, P):
    i, j = mono
    return P[0]**i * P[1]**j


def gen_C(k, pts):
    """Generator matrix of C_𝓛(pts, k∞) on E₂."""
    return matrix(F, [[ev_mono(mo, P) for P in pts] for mo in monomials(k)])


def gen_pullback(k, pts):
    """Generator matrix of C_𝓛(pts, φ*G), G = k∞.

    𝓛(φ*G) = 𝓛(k(∞) + k(T) + k(−T)) = (x − x(T))⁻ᵏ 𝓛(3k∞).
    """
    return matrix(F, [[ev_mono(mo, P) / (P[0] + xT)**k for P in pts]
                      for mo in monomials(ℓ * k)])


def scaled_gram(Gm, a):
    return Gm * matrix.diagonal(F, [c**2 for c in a]) * Gm.transpose()


def hull_dim(Gm, a):
    return Gm.nrows() - scaled_gram(Gm, a).rank()


# --------------------------------------------------------------------------
# The code φ*C, G = 11(∞), and its self-dual scaling
# --------------------------------------------------------------------------
k = 11
pullback_pts = [P + m * T for m in range(ℓ) for P in Pt]
# Check: block m of the coordinates lies over D in order: φ(P̃ᵢ ⊕ [m]T) = Pᵢ
assert all(φ(pullback_pts[m * n + i]) == D[i]
           for m in range(ℓ) for i in range(n))
Gpb = gen_pullback(k, pullback_pts)
# Check: φ*C has dimension ℓk = 33
assert Gpb.rank() == ℓ * k
shift = [((m + 1) % ℓ) * n + i for m in range(ℓ) for i in range(n)]
# Check: φ*C is quasi-cyclic: translation by T (cyclic shift of the blocks) preserves it
assert Gpb.row_space() == Gpb.matrix_from_columns(shift).row_space()
print()
print("G = %d(∞): φ*C is a QC [%d, %d, %d]₆₄ code of co-index %d"
      % (k, ℓ * n, ℓ * k, ℓ * n - ℓ * k, ℓ))
print("    length %d > q + 1 = %d, the length of the longest GRS codes"
      % (ℓ * n, q + 1))

# Residue vector r spanning C_𝓛(D, 2G)^⊥; wᵢ = √rᵢ (q even)
# Check: Σ(D) = ∞ = Σ(2G), so D − 2G is principal
assert sum(D, O2) == O2
ker = gen_C(2 * k, D).right_kernel()
# Check: C_𝓛(D, 2G)^⊥ is spanned by a single vector r
assert ker.dimension() == 1
r = ker.basis()[0]
# Check: every entry of r is nonzero, so the square roots wᵢ are nonzero
assert all(c != 0 for c in r)
w = vector(F, [c.sqrt() for c in r])
C = gen_C(k, D)
# Check: C has dimension k = 11 and w ⋆ C is self-orthogonal, hence self-dual
assert C.rank() == k and scaled_gram(C, w) == 0
w_φ = [w[i] for m in range(ℓ) for i in range(n)]
# Check: w^φ ⋆ φ*C is self-orthogonal, hence self-dual (dimension 33 = 66/2)
assert scaled_gram(Gpb, w_φ) == 0
print("Σ(D) = ∞ = Σ(2G), so D − 2G is principal")
print("w = (√rᵢ)ᵢ, where r spans C_𝓛(D, 2G)^⊥ (normalised as Sage returns it):")
print(fvec(w))
print("w ⋆ C and w^φ ⋆ φ*C are self-dual")

# --------------------------------------------------------------------------
# Hull variation (Theorem 25)
# --------------------------------------------------------------------------
λ = α                                           # λ² ≠ 1


def scaling(lams):
    """w^φ with fibre i additionally scaled by lams[i] (default 1)."""
    return [w[i] * lams.get(i, 1) for m in range(ℓ) for i in range(n)]


print()
print("λ = %s" % fe(λ))
hulls = [hull_dim(Gpb, scaling({i: λ for i in range(s)})) for s in range(k)]
# Check: Theorem 25(b): for s < k, the hull dimension is 3(k − s)
assert hulls == [ℓ * (k - s) for s in range(k)]
print("Scaling the fibres over P₁, …, Pₛ by λ, s = 0, …, %d:" % (k - 1))
print("    hull dimensions", hulls, "= 3(11 − s)")


def case(I):
    """The case of Theorem 25(c) for |I| = k, decided by Σ(D_I) ⊖ Σ(G)."""
    σ = sum((D[i] for i in I), O2)
    if σ not in ker_φ_dual:
        return "LCD"
    return "j₀ ≠ 0" if σ != O2 else "j₀ = 0"


# First subsets I (in lexicographic order) of each case
examples = {}
for I in combinations(range(n), k):
    examples.setdefault(case(I), I)
    if len(examples) == 3:
        break
# Check: all three cases of Theorem 25(c) occur among the subsets I of size k
assert len(examples) == 3

print("s = %d: the three cases of Theorem 25(c)" % k)
extra = []
for c in ("LCD", "j₀ ≠ 0", "j₀ = 0"):
    I = examples[c]
    σ = sum((D[i] for i in I), O2)
    h = hull_dim(Gpb, scaling({i: λ for i in I}))
    print("    %-8s I = %s, Σ(D_I) = %s" % (c + ":", [i + 1 for i in I], fpt(σ)))
    print("             all λᵢ = λ: hull dimension %d" % h)
    if c == "LCD":
        # Check: Theorem 25(c), LCD case: hull dimension 0
        assert h == 0
    elif c == "j₀ ≠ 0":
        # Check: Theorem 25(c), j₀ ≠ 0 case: hull dimension 2
        assert h == 2
    else:
        # Check: Theorem 25(c), j₀ = 0 case with equal λᵢ: (eqn:hull-one) holds, hull dimension 2
        assert h == 2
        # one λᵢ₀ = λ′ with 1, λ², λ′² distinct: (eqn:hull-one) fails
        λ1 = α**2
        # Check: 1, λ², λ′² are distinct, as in Corollary 26
        assert len({F(1), λ**2, λ1**2}) == 3
        lams = {i: λ for i in I}
        lams[I[0]] = λ1
        h1 = hull_dim(Gpb, scaling(lams))
        # Check: Theorem 25(c), j₀ = 0 case with one λᵢ₀ different: (eqn:hull-one) fails, hull dimension 1
        assert h1 == 1
        print("             λ_{%d} = %s, other λᵢ = λ: hull dimension %d"
              % (I[0] + 1, fe(λ1), h1))
        extra.append(h1)
    extra.append(h)

hull_set = sorted(set(hulls) | set(extra))
cs = sorted(ℓ * k - h for h in hull_set if h < ℓ * k)
print()
print("Hull dimensions h ∈ %s" % hull_set)
print("⇒ QC EAQECCs ⟦%d, c, %d; c⟧₆₄ with c = %d − h ∈ %s"
      % (ℓ * n, ℓ * n - ℓ * k, ℓ * k, cs))
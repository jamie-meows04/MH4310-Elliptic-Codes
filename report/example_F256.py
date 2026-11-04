"""
Pullback codes over 𝔽₂₅₆: Examples 18 and 27 of
"Isogeny Pullbacks of Elliptic Codes".

Coordinates of φ*D are ordered block by block: block m consists of
the points P̃ᵢ ⊕ [m]T for i = 1, …, n.
"""

from sage.all import (GF, PolynomialRing, EllipticCurve, matrix, vector,
                      floor, sqrt)

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
    """The term c·mono, e.g. α²³⁸x, xy or 1."""
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


def fvec(v, per_line=8, indent="    "):
    """A vector of field elements, several entries per line."""
    items = [fe(c) for c in v]
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
F = GF(2**8, 'α', modulus=X**8 + X**4 + X**3 + X**2 + 1)
α = F.gen()
q = F.order()


def alog(c):
    """Exponent e with c = αᵉ."""
    return c.log(α)


header("Example 18")
print("𝔽₂₅₆ = 𝔽₂(α), minimal polynomial of α:  %s = 0"
      % fpoly(F.modulus()))
print("    i.e. α⁸ = %s" % fpoly(F.modulus() - X**8))
# Check: α is a primitive element of 𝔽₂₅₆ (order 255)
assert α.multiplicative_order() == q - 1

E1 = EllipticCurve(F, [1, 0, 0, 0, 1])
pts1 = E1.points()                              # includes ∞
O1 = E1(0)
N = len(pts1)
M = q + 1 + floor(2 * sqrt(q))
# Check: N = #E₁(𝔽₂₅₆) = 288, and the Hasse–Weil bound is M = 289
assert N == 288 and M == 289
print(fcurve(E1, "E₁"))
print("N = #E₁(𝔽₂₅₆) = %d" % N)
print("M = q + 1 + ⌊2√q⌋ = %d   (Hasse–Weil bound; N = M − 1)" % M)

two_torsion = [P for P in pts1 if P != O1 and 2 * P == O1]
# Check: (0, 1) is the only rational point of order 2 on E₁
assert two_torsion == [E1(0, 1)]
H = set(2 * P for P in pts1)                    # H = [2]E₁(𝔽_q)
# Check: H = [2]E₁(𝔽₂₅₆) has index 2
assert len(H) == N // 2
print("Rational points of order 2 on E₁: %s" % fpts(two_torsion, indent=""))
print("H = [2]E₁(𝔽₂₅₆) has order %d and index %d" % (len(H), N // len(H)))

# Kernel K = ⟨T⟩ of order ℓ = 3, and Vélu's isogeny φ: E₁ → E₂
T = E1(α**238, α**176)
# Check: T has order ℓ = 3
assert T.order() == 3
ℓ = 3
K = [O1, T, 2 * T]
E2 = EllipticCurve(F, [1, 0, 0, α**238, α**51])
φ = E1.isogeny(T, codomain=E2)                  # Vélu's model is exactly E₂
O2 = E2(0)
print("T = %s, of order ℓ = %d" % (fpt(T), T.order()))
print("K = ⟨T⟩ = {%s}" % ", ".join(fpt(P) for P in K))
print("φ: E₁ → E₂ is Vélu's isogeny with kernel K, where")
print("    " + fcurve(E2, "E₂"))

# --------------------------------------------------------------------------
# Eigenfunction u, divisor A = (R) − (∞), and ker φ̂
# --------------------------------------------------------------------------
ζ = α**85
# Check: ζ is a primitive cube root of unity
assert ζ**3 == 1 and ζ != 1
print("ζ = %s, of order %d" % (fe(ζ), ζ.multiplicative_order()))

xT = T[0]
u_c1, u_c0 = α**123, α**34                      # u = (y + c₁x + c₀)/(x + x(T))


def u(P):
    """u = (y + α¹²³x + α³⁴)/(x + α²³⁸), for P ∉ K."""
    x, y = P[0], P[1]
    return (y + u_c1 * x + u_c0) / (x + xT)


affine_off_K = [P for P in pts1 if P not in K]
# Check: u is a ζ-eigenfunction of τ = τ_T*: u(P ⊕ T) = ζ·u(P)
assert all(u(P + T) == ζ * u(P) for P in affine_off_K if P + T not in K)
print("u = (y + %sx + %s)/(x + %s), with τu = ζu for τ = τ_T*"
      % (fe(u_c1), fe(u_c0), fe(xT)))

R = E2(α**51, α**20)
# Check: R has order ℓ = 3
assert R.order() == 3
# Check: the zeros of u away from K are exactly the points over R, so div u = φ*A
assert all((u(P) == 0) == (φ(P) == R) for P in affine_off_K)
fibre_R = [P for P in pts1 if φ(P) == R]
# Check: the fibre over R consists of ℓ rational points
assert len(fibre_R) == 3
print("R = %s, and div u = φ*A with A = (R) − (∞)" % fpt(R))
print("    zeros of u: φ⁻¹(R) = {%s}" % ", ".join(fpt(P) for P in fibre_R))

# E₁[3] is rational and ker φ̂ = φ(E₁[3]) = ⟨R⟩
E1_3 = [P for P in pts1 if 3 * P == O1]
# Check: E₁[3] ⊆ E₁(𝔽₂₅₆): all 9 points of order dividing 3 are rational
assert len(E1_3) == 9
ker_φ_dual = set(φ(P) for P in E1_3)
# Check: φ(E₁[3]) = ⟨R⟩
assert ker_φ_dual == {O2, R, 2 * R}
# (Sage cannot compute φ̂ in characteristic 2, so check ⟨R⟩ = ker φ̂ directly:
#  ψ: E₂ → E₂/⟨R⟩ satisfies ψ∘φ ≅ [3], i.e. ψ = φ̂ up to isomorphism.)
ψ = E2.isogeny(R)
# Check: the isogeny ψ with kernel ⟨R⟩ maps back to (a curve isomorphic to) E₁
assert ψ.codomain().is_isomorphic(E1)
# Check: ker(ψ∘φ) ∩ E₁(𝔽₂₅₆) = E₁[3], so ψ∘φ ≅ [3] and ker φ̂ = ⟨R⟩
assert [P for P in pts1 if ψ(φ(P)) == ψ.codomain()(0)] == E1_3
R2 = 2 * R
# Check: [2]R = (α⁵¹, α⁶⁵)
assert R2 == E2(α**51, α**65)
print("[2]R = %s" % fpt(R2))
print("E₁[3] ⊆ E₁(𝔽₂₅₆) (%d points), and ker φ̂ = φ(E₁[3]) = ⟨R⟩ = {%s}"
      % (len(E1_3), ", ".join(fpt(P) for P in (O2, R, R2))))

# --------------------------------------------------------------------------
# D = φ(Q ⊕ H) and the fibres over it
# --------------------------------------------------------------------------
fibres = {}
for P in pts1:
    if P not in H:                              # P ∈ Q ⊕ H
        fibres.setdefault(φ(P), []).append(P)


def point_key(P):
    return (alog(P[0]), alog(P[1]) if P[1] != 0 else -1)


Q = min((P for P in pts1 if P not in H), key=point_key)
D = sorted(fibres, key=point_key)
n = len(D)
# Check: D has n = 48 points, and each fibre over D consists of ℓ rational points
assert n == 48 and all(len(fibres[P]) == ℓ for P in D)
# Check: ∞ ∉ supp D
assert O2 not in D

D_exps = sorted(set(alog(P[0]) for P in D))
# Check: the x-coordinates of D are the αᵉ listed in the paper
assert D_exps == [5, 11, 57, 63, 79, 80, 88, 94, 109, 120, 133, 135,
                  147, 154, 169, 176, 188, 190, 203, 214, 229, 235, 243, 244]
# Check: R, [2]R ∉ supp D
assert R not in D and R2 not in D
print("Q = %s ∉ H" % fpt(Q))
print("D = φ(Q ⊕ H) has n = %d points P₁, …, P%s:" % (n, sub(n)))
print(fpts(D))
print("    x-coordinates αᵉ for e ∈ %s" % D_exps)
print("R, [2]R ∉ supp D")

# Representatives P̃ᵢ and the block ordering of φ*D
Pt = [fibres[P][0] for P in D]
pullback_pts = [P + m * T for m in range(ℓ) for P in Pt]
# Check: block m of the coordinates lies over D in order: φ(P̃ᵢ ⊕ [m]T) = Pᵢ
assert all(φ(pullback_pts[m * n + i]) == D[i]
           for m in range(ℓ) for i in range(n))
print("Representatives P̃ᵢ ∈ φ⁻¹(Pᵢ), i = 1, …, %d:" % n)
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


def h_R(P):
    """Function with simple poles at R and ∞: (y + y(−R))/(x + x(R))."""
    xR, yR = R[0], R[1]
    return (P[1] + yR + xR) / (P[0] + xR)


def gen_C(k, pts):
    """Generator matrix of C_𝓛(pts, k∞) on E₂."""
    return matrix(F, [[ev_mono(mo, P) for P in pts] for mo in monomials(k)])


def gen_B(k, j):
    """Generator matrix of Bⱼ = C_𝓛(D, Gⱼ), Gⱼ = G + jA = (k − j)(∞) + j(R).

    Basis: 𝓛((k − j)∞) together with h_R, …, h_Rʲ. As h_Rᵉ has a pole of
    order e at ∞ as well as at R, this needs k − j ≥ j.
    """
    # Check: k − j ≥ j, so h_R, …, h_Rʲ lie in 𝓛(Gⱼ)
    assert k - j >= j
    rows = [[ev_mono(mo, P) for P in D] for mo in monomials(k - j)]
    rows += [[h_R(P)**e for P in D] for e in range(1, j + 1)]
    Bj = matrix(F, rows)
    # Check: Bⱼ has dimension k
    assert Bj.rank() == k
    return Bj


def gen_pullback(k, pts):
    """Generator matrix of C_𝓛(pts, φ*G), G = k∞.

    𝓛(φ*G) = 𝓛(k(∞) + k(T) + k(−T)) = (x − x(T))⁻ᵏ 𝓛(3k∞).
    """
    return matrix(F, [[ev_mono(mo, P) / (P[0] + xT)**k for P in pts]
                      for mo in monomials(ℓ * k)])


# --------------------------------------------------------------------------
# Example 18: [144, 3k, 145 − 3k] QC MDS codes and matrix-product structure
# --------------------------------------------------------------------------

v = [[u(P)**j for P in Pt] for j in range(ℓ)]  # vⱼ = (u(P̃ᵢ)ʲ)ᵢ
for j in range(ℓ):
    print("v%s = (u(P̃ᵢ)%s)ᵢ:" % (sub(j), sup(j)))
    print(fvec(v[j]))

V = matrix(F, [[ζ**(j * m) for m in range(ℓ)] for j in range(ℓ)])
print("V = (ζʲᵐ)ⱼ,ₘ₌₀² =")
for row in V.rows():
    print("    [" + ", ".join("%4s" % fe(c) for c in row) + " ]")

def matrix_product(k):
    """Generator matrix of (C₀, C₁, C₂)·V with Cⱼ = vⱼ ⋆ Bⱼ."""
    rows = []
    for j in range(ℓ):
        for r in gen_B(k, j).rows():
            cw = [r[i] * v[j][i] for i in range(n)]
            rows.append([ζ**(j * m) * cw[i]
                         for m in range(ℓ) for i in range(n)])
    return matrix(F, rows)

print()
print("G = k(∞), k odd: φ*C is a QC [%d, 3k, %d − 3k]₂₅₆ MDS code of co-index 3"
      % (ℓ * n, ℓ * n + 1))
print("Length %d = ⌊M/2⌋; Munuera's bound is (q+1)/2 + √q + κ = %.1f + κ"
      % (ℓ * n, float((q + 1) / 2 + sqrt(q))))
print("Twisted divisors Gⱼ = G + jA = (k − j)(∞) + j(R),    j = 0, 1, 2")
print("Twisted codes Cⱼ = vⱼ ⋆ C_𝓛(D, Gⱼ),   "
      "[%d, k, %d − k]₂₅₆ MDS codes" % (n, n + 1))

def check_example_18(k):
    # Check: k is odd (Construction 17(a)), at least 5 (for gen_B), and less than n
    assert k % 2 == 1 and 5 <= k < n
    Gpb = gen_pullback(k, pullback_pts)
    # Check: φ*C has dimension ℓk
    assert Gpb.rank() == ℓ * k

    # translation by T permutes coordinates block-cyclically
    # test quasi-cyclic structure of φ*C
    shift = [((m + 1) % ℓ) * n + i for m in range(ℓ) for i in range(n)]
    Gshift = Gpb.matrix_from_columns(shift)
    # Check: φ*C is quasi-cyclic: the cyclic shift of the blocks preserves its row space
    assert Gpb.row_space() == Gshift.row_space()

    # test exact matrix-product decomposition
    MP = matrix_product(k)
    # Check: φ*C = (C₀, C₁, C₂)·V exactly (equal row spaces)
    assert MP.row_space() == Gpb.row_space()
    
    print("k = %2d: φ*C is a QC [%d, %d] code of co-index 3, "
          "equal to (C₀, C₁, C₂)·V" % (k, ℓ * n, ℓ * k))


print("\nChecks:")
for k in (5, 13, 23):
    check_example_18(k)

# MDS check for the smallest case k = 1 ([144, 3, 142]):
# every 3 columns of the generator matrix are independent.
# (Slow-ish: about 4.8 × 10⁵ determinants. Set RUN_MDS_CHECK = True to run.)
RUN_MDS_CHECK = False
if RUN_MDS_CHECK:
    from itertools import combinations
    G1 = gen_pullback(1, pullback_pts)
    cols = G1.columns()
    # Check: every 3 columns of the k = 1 generator matrix are independent, so the code is MDS
    assert all(matrix([cols[a], cols[b], cols[c]]).det() != 0
               for a, b, c in combinations(range(ℓ * n), 3))
    print("k = 1: [144, 3, 142] MDS confirmed")

# The 95 affine points of φ(E₁(𝔽_q)) give QC [285, 3k] codes
image = set(φ(P) for P in pts1)
# Check: φ(E₁(𝔽₂₅₆)) has 96 points (95 affine points and ∞)
assert len(image) == 96 and O2 in image
print()
print("φ(E₁(𝔽₂₅₆)) has %d affine points: QC [%d, 3k]₂₅₆ codes for 1 ≤ k ≤ %d"
      % (len(image) - 1, ℓ * (len(image) - 1), len(image) - 2))

# --------------------------------------------------------------------------
# Example 27: self-dual QC [138, 69, 70] MDS code and its hull variation
# --------------------------------------------------------------------------
header("Example 27")
k = 23
removed = [P for P in D if alog(P[0]) == 5]
Dp = [P for P in D if alog(P[0]) != 5]          # remove the two points x = α⁵
n2 = len(Dp)
# Check: D′ has 46 points
assert n2 == 46
# Check: Σ(D′) = ∞ = Σ(2G), so D′ − 2G is principal
assert sum(Dp, O2) == O2
print("k = %d, G = %d(∞)" % (k, k))
print(r"D′ = D \ {%s}, with %d points" % (", ".join(fpt(P) for P in removed),
                                             n2))
print("Σ(D′) = %s = Σ(2G), so D′ − 2G is principal" % fpt(sum(Dp, O2)))
print("φ*C = C_𝓛(φ*D′, φ*G) is a QC [%d, %d, %d]₂₅₆ MDS code of co-index 3"
      % (ℓ * n2, ℓ * k, ℓ * n2 - ℓ * k + 1))

# Residue vector r spanning C_𝓛(D′, 2G)^⊥; wᵢ = √rᵢ (q even)
ker = gen_C(2 * k, Dp).right_kernel()
# Check: C_𝓛(D′, 2G)^⊥ is spanned by a single vector r
assert ker.dimension() == 1
r = ker.basis()[0]
# Check: every entry of r is nonzero, so the square roots wᵢ are nonzero
assert all(c != 0 for c in r)
w = vector(F, [c.sqrt() for c in r])

C = gen_C(k, Dp)
# Check: C has dimension k = 23
assert C.rank() == k


def scaled_gram(Gm, a):
    return Gm * matrix.diagonal(F, [c**2 for c in a]) * Gm.transpose()


# Check: w ⋆ C is self-orthogonal, hence self-dual (dimension 23 = 46/2)
assert scaled_gram(C, w) == 0

# Pullback along D′, coordinates in blocks as before
Pt2 = [fibres[P][0] for P in Dp]
pts2 = [P + m * T for m in range(ℓ) for P in Pt2]
Gpb = gen_pullback(k, pts2)
# Check: φ*C has dimension ℓk = 69
assert Gpb.rank() == ℓ * k
w_φ = [w[i] for m in range(ℓ) for i in range(n2)]
# Check: w^φ ⋆ φ*C is self-orthogonal, hence self-dual (dimension 69 = 138/2)
assert scaled_gram(Gpb, w_φ) == 0
print("w = (√rᵢ)ᵢ, where r spans C_𝓛(D′, 2G)^⊥ (normalised as Sage returns it):")
print(fvec(w))
print("w ⋆ C and w^φ ⋆ φ*C are self-dual")


def hull_dim(Gm, a):
    return Gm.nrows() - scaled_gram(Gm, a).rank()


λ = α                                           # λ² ≠ 1
hulls = []
for s in range(k + 1):
    I = set(range(s))                           # scale the first s fibres
    a = [w[i] * (λ if i in I else 1) for m in range(ℓ) for i in range(n2)]
    hulls.append(hull_dim(Gpb, a))
# Check: Theorem 25(b): for s < k, the hull dimension is 3(k − s)
assert hulls[:k] == [ℓ * (k - s) for s in range(k)]
# Theorem 25(c): for s = k, the case is decided by whether Σ(D_I) ∈ ker φ̂
σ_I = sum(Dp[:k], O2)
if σ_I not in ker_φ_dual:
    expected = 0                                # LCD
elif σ_I != O2:
    expected = 2
else:
    expected = None                             # 1 or 2, by (eqn:hull-one)
# Check: Theorem 25(c): for s = k, the hull dimension matches the case given by Σ(D_I)
assert expected is None or hulls[k] == expected
print("λ = %s; scaling the fibres over P′₁, …, P′ₛ by λ" % fe(λ))
print("Hull dimensions for s = 0, …, %d:" % k, hulls)
print("    for s = %d: Σ(D_I) = %s" % (k, fpt(σ_I)))
print("⇒ EAQECCs ⟦%d, 3s, %d; 3s⟧₂₅₆ for 1 ≤ s ≤ %d"
      % (ℓ * n2, ℓ * n2 - ℓ * k + 1, k))
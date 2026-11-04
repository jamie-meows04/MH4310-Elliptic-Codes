"""
Pullback codes over 𝔽₂₅₁: Example 21 of MH4310 report
"Elliptic Curve Codes: Recent Constructions and Isogeny-Induced
Transformations".

Coordinates of φ*D are ordered block by block: block m consists of
the points P̃ᵢ ⊕ [m]T for i = 1, …, n.
"""

from sage.all import GF, EllipticCurve, matrix, floor, sqrt

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
    """A field element as an integer in {0, …, q − 1}."""
    return str(int(c))


def fpt(P):
    """A point as (x, y), or ∞."""
    if P.is_zero():
        return "∞"
    return "(%s, %s)" % (fe(P[0]), fe(P[1]))


def fterm(c, mono):
    """The term c·mono, e.g. 231x, x or 13."""
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


def fvec(v, per_line=14, indent="    "):
    """A vector of field elements, several entries per line."""
    items = ["%3s" % fe(c) for c in v]
    lines = [", ".join(items[i:i + per_line])
             for i in range(0, len(items), per_line)]
    return indent + ("\n" + indent).join(lines)


def fpts(pts, per_line=7, indent="    "):
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


def point_key(P):
    """Lexicographic order on affine points, by integer coordinates."""
    return (int(P[0]), int(P[1]))


# --------------------------------------------------------------------------
# Field and curves
# --------------------------------------------------------------------------
q = 251
F = GF(q)

header("Example 18 over 𝔽₂₅₁")
print("𝔽₂₅₁ = ℤ/251ℤ")

E1 = EllipticCurve(F, [1, 100])                 # y² = x³ + x + 100
pts1 = E1.points()                              # includes ∞
O1 = E1(0)
N = len(pts1)
M = q + 1 + floor(2 * sqrt(q))
# Check: N = #E₁(𝔽₂₅₁) = 280, and the Hasse–Weil bound is M = 283
assert N == 280 and M == 283
print(fcurve(E1, "E₁"))
print("N = #E₁(𝔽₂₅₁) = %d" % N)
print("M = q + 1 + ⌊2√q⌋ = %d   (Hasse–Weil bound; N = M − 3)" % M)

two_torsion = [P for P in pts1 if P != O1 and 2 * P == O1]
# Check: (199, 0) is the only rational point of order 2 on E₁
assert two_torsion == [E1(199, 0)]
H = set(2 * P for P in pts1)                    # H = [2]E₁(𝔽_q)
# Check: H = [2]E₁(𝔽₂₅₁) has index 2
assert len(H) == N // 2
print("Rational points of order 2 on E₁: %s" % fpts(two_torsion, indent=""))
print("H = [2]E₁(𝔽₂₅₁) has order %d and index %d" % (len(H), N // len(H)))

# Kernel K = ⟨T⟩ of order ℓ = 5, and Vélu's isogeny φ: E₁ → E₂
T = E1(65, 14)
# Check: T has order ℓ = 5
assert T.order() == 5
ℓ = 5
K = [m * T for m in range(ℓ)]
E2 = EllipticCurve(F, [231, 13])                # y² = x³ + 231x + 13
O2 = E2(0)
# Sage may realise φ only up to an automorphism of E₂ (here ±1). Fix the
# sign so that φ is Vélu's normalised isogeny, as in the paper: φ(0, 10) =
# (71, 109) by Vélu's formulae.
φ_sage = E1.isogeny(T, codomain=E2)
sign = 1 if φ_sage(E1(0, 10)) == E2(71, 109) else -1


def φ(P):
    return sign * φ_sage(P)


# Check: φ is Vélu's normalised isogeny (not its negative)
assert φ(E1(0, 10)) == E2(71, 109)
print("T = %s, of order ℓ = %d" % (fpt(T), T.order()))
print("K = ⟨T⟩ = {%s}" % ", ".join(fpt(P) for P in K))
print("φ: E₁ → E₂ is Vélu's isogeny with kernel K, where")
print("    " + fcurve(E2, "E₂"))

# --------------------------------------------------------------------------
# Eigenfunction u, divisor A = (R) − (∞), and ker φ̂
# --------------------------------------------------------------------------
ζ = F(20)
# Check: ℓ divides q − 1, so 𝔽₂₅₁ contains a primitive ℓ-th root of unity
assert (q - 1) % ℓ == 0
# Check: ζ = 20 is a primitive ℓ-th root of unity
assert ζ.multiplicative_order() == ℓ
print("ζ = %s, of order %d" % (fe(ζ), ζ.multiplicative_order()))

# u = (xy + 195x² + 99y + 215x + 232)/((x − 65)(x − 204))
u_num = {"xy": F(1), "x²": F(195), "y": F(99), "x": F(215), "": F(232)}
u_den_roots = (F(65), F(204))                   # x(T) and x([2]T)
# Check: the roots of the denominator of u are the x-coordinates of the affine points of K
assert {P[0] for P in K[1:]} == set(u_den_roots)


def u_numerator(P):
    x, y = P[0], P[1]
    return (u_num["xy"] * x * y + u_num["x²"] * x**2 + u_num["y"] * y
            + u_num["x"] * x + u_num[""])


def u(P):
    """u = (xy + 195x² + 99y + 215x + 232)/((x − 65)(x − 204)), for P ∉ K."""
    x = P[0]
    return u_numerator(P) / ((x - u_den_roots[0]) * (x - u_den_roots[1]))


affine_off_K = [P for P in pts1 if P not in K]
# Check: u is a ζ-eigenfunction of τ = τ_T*: u(P ⊕ T) = ζ·u(P)
assert all(u(P + T) == ζ * u(P) for P in affine_off_K if P + T not in K)
num_str = " + ".join(t for t in (fterm(u_num[m], m) for m in u_num) if t)
print("u = (%s)/((x − %s)(x − %s)), with τu = ζu for τ = τ_T*"
      % (num_str, fe(u_den_roots[0]), fe(u_den_roots[1])))

# Poles: the denominator vanishes at the four affine points of K, and the
# numerator does not, so u has simple poles there; u also has a simple
# pole at ∞ (numerator of order 5, denominator of order 4).
# Check: the numerator of u does not vanish on K ∖ {∞}, so u has simple poles at the points of K
assert all(u_numerator(P) != 0 for P in K[1:])
print("    poles of u: the five points of K")

R = E2(141, 209)
# Check: R has order ℓ
assert R.order() == ℓ
# Check: the rational zeros of u away from K are exactly the rational points over R
assert all((u(P) == 0) == (φ(P) == R) for P in affine_off_K)
fibre_R = [P for P in pts1 if φ(P) == R]
print("R = %s, and div u = φ*A with A = (R) − (∞)" % fpt(R))
if len(fibre_R) == ℓ:
    print("    zeros of u: φ⁻¹(R) = {%s}" % ", ".join(fpt(P) for P in fibre_R))
else:
    # Check: no point over R is rational (the fibre is defined over an extension of 𝔽₂₅₁)
    assert fibre_R == []
    print("    zeros of u: φ⁻¹(R), which has no 𝔽₂₅₁-rational points")

# ker φ̂ = ⟨R⟩ (φ̂ is available in odd characteristic ≠ 3)
φ_dual = φ_sage.dual()                        # same kernel as for φ
# Check: φ̂ maps E₂ back to (a curve isomorphic to) E₁
assert φ_dual.codomain().is_isomorphic(E1)
ker_φ_dual = [m * R for m in range(ℓ)]
# Check: ⟨R⟩ ⊆ ker φ̂; both have order ℓ, so ker φ̂ = ⟨R⟩
assert all(φ_dual(P) == φ_dual.codomain()(0) for P in ker_φ_dual)
print("ker φ̂ = ⟨R⟩ = {%s}" % ", ".join(fpt(P) for P in ker_φ_dual))

# --------------------------------------------------------------------------
# D: the points of E₂ whose fibres lie in E₁(𝔽_q) ∖ H
# --------------------------------------------------------------------------
fibres = {}
for P in pts1:
    if P not in H:                              # P ∈ Q ⊕ H
        fibres.setdefault(φ(P), []).append(P)

D_list = [(5, 17), (23, 26), (56, 8), (60, 31), (77, 52), (104, 24),
          (142, 124), (143, 67), (161, 100), (222, 7), (224, 85),
          (225, 5), (228, 75), (232, 78)]
D = []
for (x, y) in D_list:                           # (x, y) before (x, −y)
    D += [E2(x, y), E2(x, -y)]
n = len(D)
# Check: D has 28 points, and they are exactly the images φ(P) of P ∈ E₁(𝔽₂₅₁) ∖ H
assert n == 28 and set(D) == set(fibres)
# Check: each fibre over D consists of ℓ rational points
assert all(len(fibres[P]) == ℓ for P in D)
# Check: ∞ ∉ supp D and ⟨R⟩ ∩ supp D = ∅
assert O2 not in D and not any(P in D for P in ker_φ_dual)
print("D has n = %d points P₁, …, P%s, ordered with (x, y) before (x, −y):"
      % (n, sub(n)))
print(fpts(D))
print("⟨R⟩ ∩ supp D = ∅")

# Representatives P̃ᵢ: the lexicographically smallest point of φ⁻¹(Pᵢ)
Pt = [min(fibres[P], key=point_key) for P in D]
pullback_pts = [P + m * T for m in range(ℓ) for P in Pt]
# Check: block m of the coordinates lies over D in order: φ(P̃ᵢ ⊕ [m]T) = Pᵢ
assert all(φ(pullback_pts[m * n + i]) == D[i]
           for m in range(ℓ) for i in range(n))
print("Representatives P̃ᵢ ∈ φ⁻¹(Pᵢ), lexicographically smallest:")
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
    """Function with simple poles at R and ∞: (y + y(R))/(x − x(R))."""
    xR, yR = R[0], R[1]
    return (P[1] + yR) / (P[0] - xR)


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

    With d = (x − 65)(x − 204), div d = Σ_{T′ ∈ K∖{∞}} (T′) − 4(∞), so
    𝓛(φ*G) = 𝓛(k Σ_{T′ ∈ K} (T′)) = d⁻ᵏ 𝓛(5k∞).
    """
    def d(P):
        return (P[0] - u_den_roots[0]) * (P[0] - u_den_roots[1])
    return matrix(F, [[ev_mono(mo, P) / d(P)**k for P in pts]
                      for mo in monomials(ℓ * k)])


# --------------------------------------------------------------------------
# Parameters and matrix-product structure
# --------------------------------------------------------------------------

v = [[u(P)**j for P in Pt] for j in range(ℓ)]  # vⱼ = (u(P̃ᵢ)ʲ)ᵢ
for j in range(ℓ):
    print("v%s = (u(P̃ᵢ)%s)ᵢ:" % (sub(j), sup(j)))
    print(fvec(v[j]))

V = matrix(F, [[ζ**(j * m) for m in range(ℓ)] for j in range(ℓ)])
V_paper = matrix(F, [[1, 1, 1, 1, 1],
                     [1, 20, 149, 219, 113],
                     [1, 149, 113, 20, 219],
                     [1, 219, 20, 113, 149],
                     [1, 113, 219, 149, 20]])
# Check: V matches the matrix in the paper
assert V == V_paper
print("V = (ζʲᵐ)ⱼ,ₘ₌₀⁴ =")
for row in V.rows():
    print("    [" + ", ".join("%3s" % fe(c) for c in row) + " ]")

print()
print("G = k(∞), k odd: φ*C is a QC [%d, %dk, %d − %dk]₂₅₁ MDS code "
      "of co-index %d" % (ℓ * n, ℓ, ℓ * n + 1, ℓ, ℓ))
print("Twisted divisors Gⱼ = G + jA = (k − j)(∞) + j(R), j = 0, …, %d" % (ℓ - 1))
print("Twisted codes Cⱼ = vⱼ ⋆ C_𝓛(D, Gⱼ),   "
      "[%d, k, %d − k]₂₅₁ MDS codes" % (n, n + 1))

def matrix_product(k):
    """Generator matrix of (C₀, …, C₄)·V with Cⱼ = vⱼ ⋆ Bⱼ."""
    rows = []
    for j in range(ℓ):
        for r in gen_B(k, j).rows():
            cw = [r[i] * v[j][i] for i in range(n)]
            rows.append([ζ**(j * m) * cw[i]
                         for m in range(ℓ) for i in range(n)])
    return matrix(F, rows)


def check(k):
    # Check: k is odd (Construction 17(a)), at least 2(ℓ − 1) (for gen_B), and less than n
    assert k % 2 == 1 and 2 * (ℓ - 1) <= k < n
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
    # Check: φ*C = (C₀, …, C₄)·V exactly (equal row spaces)
    assert MP.row_space() == Gpb.row_space()

    print("k = %2d: φ*C is a QC [%d, %d, %d] code of co-index %d, "
          "equal to (C₀, …, C₄)·V"
          % (k, ℓ * n, ℓ * k, ℓ * n - ℓ * k + 1, ℓ))


print("\nChecks:")
for k in (9, 15, 25):
    check(k)

# The 55 affine points of φ(E₁(𝔽_q)) give QC [275, 5k] codes
image = set(φ(P) for P in pts1)
# Check: φ(E₁(𝔽₂₅₁)) has 56 points (55 affine points and ∞)
assert len(image) == 56 and O2 in image
print("φ(E₁(𝔽₂₅₁)) has %d affine points: QC [%d, %dk]₂₅₁ codes for 1 ≤ k ≤ %d"
      % (len(image) - 1, ℓ * (len(image) - 1), ℓ, len(image) - 2))
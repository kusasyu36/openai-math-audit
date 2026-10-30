"""Search for a rational point on E_l: y^2=x(x-l)(x+3l) in Kummer class (beta,1)
via the 2-covering:  x = beta*Z^2/W^2,  x - l = Y^2/W^2,  x + 3l = beta*U^2/W^2.
Parametrize the conic beta Z^2 - Y^2 = l W^2 from one small point, then search the line parameter.
"""
import sys
from math import isqrt, gcd
from fractions import Fraction

def is_sq(n):
    return n >= 0 and isqrt(n)**2 == n

def conic_point(beta, l, B=3000):
    for W in range(1, B):
        lw = l*W*W
        # beta Z^2 = Y^2 + l W^2 ; Z >= sqrt(l/beta) W
        z0 = isqrt(lw//beta)
        for Z in range(max(1, z0), z0 + B):
            v = beta*Z*Z - lw
            if v < 0: continue
            if is_sq(v):
                return (Z, isqrt(v), W)
    return None

def search(beta, l, B=400):
    pt = conic_point(beta, l)
    if pt is None: return None
    Z0, Y0, W0 = pt
    # lines through (Z0:Y0:W0): (Z,Y,W) = t*(Z0,Y0,W0) + (a,b,c) direction; use standard param:
    # points (Z0 + p*u, Y0 + q*u, W0 + r*u)... simpler: parametrize with slope (m,n): Z = Z0*k + m, etc.
    # Use generic formula: for direction D=(dz,dy,dw), second intersection: s = -2 B(P,D)/Q(D)
    def Qf(z, y, w): return beta*z*z - y*y - l*w*w
    def Bf(P, D): return beta*P[0]*D[0] - P[1]*D[1] - l*P[2]*D[2]
    P = (Z0, Y0, W0)
    for a in range(-B, B+1):
        for b in range(0, B+1):
            for D in ((a, b, 1), (a, 1, b)):
                qd = Qf(*D)
                if qd == 0: continue
                bd = Bf(P, D)
                # point = qd*P - 2*bd*D  (scaled)
                Z = qd*P[0] - 2*bd*D[0]; Y = qd*P[1] - 2*bd*D[1]; W = qd*P[2] - 2*bd*D[2]
                if W == 0 or Z == 0: continue
                g = gcd(gcd(abs(Z), abs(Y)), abs(W)); Z//=g; Y//=g; W//=g
                val = beta*(beta*Z*Z + 3*l*W*W)
                if is_sq(val):
                    x = Fraction(beta*Z*Z, W*W)
                    f = x*(x-l)*(x+3*l)
                    assert is_sq(f.numerator) and is_sq(f.denominator)
                    return x
    return None

cases = [(171011, 3977), (1751, 17), (5963, 89), (14663, 473), (31691, 473), (75107, 3953),
         (1271, 41), (23579, 1241), (323, 17), (1763, 41), (41819, 2201), (47615, 89),
         (47027, 41), (161891, 89), (226703, 2201), (169847, 1649), (766379, 89), (411983, 737),
         (397271, 56753), (2111603, 1577)]
found = 0
for l, beta in cases:
    x = search(beta, l, B=int(sys.argv[1]) if len(sys.argv) > 1 else 60)
    print(l, beta, "x =", x)
    found += x is not None
print("found", found, "of", len(cases))

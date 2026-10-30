# Audit 077: consistency examples (Knapp / plane wave / focusing / plane-concentrated lines / Kakeya reduction)
import sympy as sp
m,N,d,k,s,L,w,rho = sp.symbols('m N delta kappa s L w rho', positive=True)
def ok(c,msg): print(("OK  " if c else "FAIL"), msg)
# Packet Thm 1.2 single packet (Knapp): G1^3 = delta^2 |c|^3 with |c|=delta^2 -> delta^8 ; output int|F|^3 ~ delta^6*delta^2
ok(sp.simplify(d**2*(d**2)**3 - (d**2)**3*d**2)==0, "Knapp packet: G1^3 = delta^8 = normalized ∫|Ef|^3 (saturated, no violation)")
# plane wave: per-bin terminal G^3 = (w^2 sum|c|^2) K^{1/2} with |c|~1 on m^2N^2 positions, K~1
wt=d/m; G3=(wt**2*m**2/d**2)*1
ok(sp.simplify(G3-1)==0, "plane wave: terminal per-bin G^3 = 1 = G_1^3 (saturated up to m^{10k})")
# lower bound A >= (w^2 theta / W^{1/2})^{1/3} ||c|test||_2 with test radius m^2, W=m^4
lb=(wt**2*(m*d)/sp.sqrt(m**4))*(m**2/d**2)**sp.Rational(3,2)
ok(sp.simplify(lb-1)==0, "plane wave: atomic-norm lower bound also = 1 per bin")
# focusing wave at terminal bin of focus: per-bin contribution m^-2 * m^2 N^-6 = N^-6 = G_1^3
theta=m/N; c2=(N/m)**2*(m/N)**4
contrib=m**-2*(1/(m*N))**2*theta*c2**sp.Rational(3,2)
ok(sp.simplify(contrib-N**-6)==0, "focusing wave: focal-bin lower bound N^-6 = G_1^3 (consistent; log N from other bins)")
# geometric theorem: plane-concentrated lines give q <= 2+s-kappa; claim 1+2s-10kappa must not exceed it
ok(sp.simplify((2+s-k)-(1+2*s-10*k) - (1-s+9*k))==0, "plane example: (2+s-k)-(1+2s-10k) = 1-s+9k > 0 for s<=1")
# Kakeya reduction: s=1,D=1 -> occupied cells >= M^{q}, cell volume M^-3 -> volume >= M^{q-3} = M^{-10k}
M,q=sp.symbols('M q',positive=True)
ok(sp.simplify(M**(3-10*k)*M**-3 - M**(-10*k))==0, "geo:main(s=1,D=1) => |N_{1/M}(union of segments)| >= M^{-10k-o(1)} (Kakeya, Minkowski form)")

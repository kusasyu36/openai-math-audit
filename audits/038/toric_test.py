# Continuum toric sanity test of Section 5 (Prop 5.1).
# P = conv{(0,0),(K,0),(1,1),(0,1)} (smooth polygon, Hirzebruch surface), x = vertex (0,0).
# Monomial tests v=(a,b) on toric coords (u,w); A(v)=a+b; Phi(a,b)=exp((a+b)/t)*mean_P exp(-a u - b w).
import mpmath as mp
mp.mp.dps=20
def lam(K,u): return 1 if u<=1 else (K-u)/(K-1)
def mean(K,a,b):
    def inner(u):
        l=lam(K,u)
        return l if b==0 else (1-mp.e**(-b*l))/b
    I=mp.quad(lambda u: mp.e**(-a*u)*inner(u),[0,1,K])
    return I/((K+1)/mp.mpf(2))
def Phi(K,t,a,b): return mp.e**((a+b)/t)*mean(K,a,b)
def golden(f,lo,hi,it=60):
    g=(mp.sqrt(5)-1)/2
    c=hi-g*(hi-lo); d=lo+g*(hi-lo)
    for _ in range(it):
        if f(c)<f(d): hi=d
        else: lo=c
        c=hi-g*(hi-lo); d=lo+g*(hi-lo)
    return (lo+hi)/2
for K in [2,5,20]:
  for t in [mp.mpf('2.9'),mp.mpf('2.99')]:
    a,b=mp.mpf(1),mp.mpf(1)
    for _ in range(8):
        a=golden(lambda x: Phi(K,t,x,b),0,30)
        b=golden(lambda y: Phi(K,t,a,y),0,30)
    a0=golden(lambda x: Phi(K,t,x,0),0,30)
    Ew = mp.quad(lambda u: mp.e**(-a0*u)*lam(K,u)**2/2,[0,1,K])/mp.quad(lambda u: mp.e**(-a0*u)*lam(K,u),[0,1,K])
    print(K,float(t),"min(a,b)=",mp.nstr(a,6),mp.nstr(b,6),"F=",mp.nstr(Phi(K,t,a,b),8),
          "| b=0 line: a0=",mp.nstr(a0,6),"F=",mp.nstr(Phi(K,t,a0,0),8),"E[w]=",mp.nstr(Ew,6),"1/t=",mp.nstr(1/t,6))

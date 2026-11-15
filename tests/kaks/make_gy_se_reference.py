"""Make tests/kaks/golden/gy_cases.se_reference.tsv: standard errors of dN and dS for method GY, computed
independently of COAL's C++ code (needs numpy; the test only reads the stored file).

    python tests/kaks/make_gy_se_reference.py build/coal-cpp

For every pair of tests/kaks/gy_cases.axt (standard genetic code) for which COAL gives both standard errors, the
log-likelihood of GY94 / F3x4 is rebuilt here from the model's definition (numpy's eigh), its Hessian in
(t, log kappa, log omega) at COAL's estimates is taken by five-point second differences of that likelihood, and the
standard errors of dS and dN follow by the delta method. COAL's own route is different (central differences of its
analytic gradient), so agreement to four digits checks both. (codeml's getSE = 1 uses a cruder Hessian and differs
from both by up to about 15 %: tests/kaks/golden/gy_cases.codeml.tsv has its values.)
"""
import csv
import io
import math
import os
import subprocess
import sys

import numpy as np

BASES="TCAG"; CODONS=[a+b+c for a in BASES for b in BASES for c in BASES]
AA="FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
def setup(a,b):
    cod=[(a[i:i+3],b[i:i+3]) for i in range(0,len(a),3)]
    f=np.zeros((3,4))
    for x,y in cod:
        for p in range(3):
            f[p][BASES.index(x[p])]+=1; f[p][BASES.index(y[p])]+=1
    f/=2*len(cod)
    sense=[c for c in range(64) if AA[c]!="*"]
    pi=np.array([f[0][c//16]*f[1][(c//4)%4]*f[2][c%4] for c in sense]); pi/=pi.sum()
    n=len(sense); ts=np.zeros((n,n)); non=np.zeros((n,n)); adj=np.zeros((n,n))
    for i,ci in enumerate(sense):
        for j,cj in enumerate(sense):
            d=[(CODONS[ci][p],CODONS[cj][p]) for p in range(3) if CODONS[ci][p]!=CODONS[cj][p]]
            if len(d)!=1: continue
            adj[i,j]=1
            if set(d[0]) in ({"T","C"},{"A","G"}): ts[i,j]=1
            if AA[ci]!=AA[cj]: non[i,j]=1
    idx={CODONS[c]:i for i,c in enumerate(sense)}
    cnt=np.zeros((n,n))
    for x,y in cod: cnt[idx[x],idx[y]]+=1
    return pi,adj,ts,non,cnt,len(cod)
def parts(M, x):
    pi,adj,ts,non,cnt,L=M; t,lk,lw=x; k,w=math.exp(lk),math.exp(lw)
    S=adj*np.where(ts>0,k,1.0)*np.where(non>0,w,1.0)
    Q=S*pi[None,:]; np.fill_diagonal(Q,-Q.sum(1)); mr=-(pi*np.diag(Q)).sum()
    return Q/mr
def lnl(M,x):
    pi,adj,ts,non,cnt,L=M; Q=parts(M,x)
    sq=np.sqrt(pi); A=sq[:,None]*Q/sq[None,:]; A=(A+A.T)/2
    d,V=np.linalg.eigh(A); P=(V*np.exp(d*x[0]))@V.T; P=P*sq[None,:]/sq[:,None]
    return (cnt*np.log(np.maximum(pi[:,None]*P,1e-300)))[cnt>0].sum()
def dsdn(M,x):
    pi,adj,ts,non,cnt,L=M; t,lk,lw=x; k,w=math.exp(lk),math.exp(lw)
    S=adj*np.where(ts>0,k,1.0)*pi[None,:]*pi[:,None]
    syn=(S*(1-non)).sum(); non0=(S*non).sum()
    dS=t*(syn+non0)/(3*(syn+w*non0)); return np.array([dS,w*dS])
def se(M,x0):
    x0=np.array(x0,float); h=np.array([2e-3*x0[0],2e-3,2e-3]); H=np.zeros((3,3))
    f0=lnl(M,x0)
    for i in range(3):
        for j in range(i,3):
            if i==j:
                e=np.zeros(3); e[i]=h[i]
                H[i,i]=(-lnl(M,x0+2*e)+16*lnl(M,x0+e)-30*f0+16*lnl(M,x0-e)-lnl(M,x0-2*e))/(12*h[i]**2)
            else:
                ei=np.zeros(3); ei[i]=h[i]; ej=np.zeros(3); ej[j]=h[j]
                H[i,j]=H[j,i]=(lnl(M,x0+ei+ej)-lnl(M,x0+ei-ej)-lnl(M,x0-ei+ej)+lnl(M,x0-ei-ej))/(4*h[i]*h[j])
    V=np.linalg.inv(-H); J=np.zeros((2,3))
    for j in range(3):
        e=np.zeros(3); e[j]=h[j]; J[:,j]=(dsdn(M,x0+e)-dsdn(M,x0-e))/(2*h[j])
    C=J@V@J.T
    return f0, math.sqrt(C[0,0]), math.sqrt(C[1,1])


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    axt = os.path.join(here, "gy_cases.axt")
    lines = open(axt).read().split("\n")
    pairs = {lines[k].split("&")[0]: (lines[k + 1], lines[k + 2]) for k in range(0, len(lines) - 3, 4)}
    out = subprocess.run([sys.argv[1], "dnds", axt, "-m", "GY", "-q", "--ml-columns", "--digits", "12"],
                         capture_output=True, text=True, check=True).stdout
    with open(os.path.join(here, "golden", "gy_cases.se_reference.tsv"), "w", newline="\n") as fh:
        fh.write("name\tlnL\tSE_dN\tSE_dS\n")
        for r in csv.DictReader(io.StringIO(out), delimiter="\t"):
            name = r["seq1"]
            if name.endswith("_mt") or "NA" in (r["SE_dN"], r["SE_dS"]):
                continue
            M = setup(*pairs[name])
            f0, ss, sn = se(M, [float(r["t"]), math.log(float(r["kappa"])), math.log(float(r["dNdS"]))])
            fh.write(f"{name}\t{f0:.6f}\t{sn:.6g}\t{ss:.6g}\n")
            print(name, f0, r["lnL"], sn, r["SE_dN"], ss, r["SE_dS"])


if __name__ == "__main__":
    main()

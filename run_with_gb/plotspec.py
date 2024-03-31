#python seg.py AET*.dat
import matplotlib
matplotlib.use("Agg")
from optparse import OptionParser
import subprocess
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import glob

def load_text(file):
        data = np.loadtxt(file, delimiter=' ')
        ft = np.array(data[:,0])
        A = np.array(data[:,1])
        E = np.array(data[:,2])
        T = np.array(data[:,3])
        return ft, A, E, T

f1, s1n, s1m, s1nm = load_text("./specfit_0_0.dat")
f2, Ar, Ai, Er = load_text("./AET_seg0_f_complex.dat")
fig1=plt.figure()
plt.plot(f1, s1nm, color="red", alpha=0.5, label="0 specfit col3")
ps=(Ar**2)+(Ai**2)
plt.plot(f2, ps, color="blue", alpha=0.5, label="ps")
#plt.ylim(-0.5e-18,0.5e-18)
plt.yscale('log')
plt.legend()
plt.xscale('log')
plt.ylim(10e-49,10e-37)
plt.xlim(10e-06,10e-01)
fig1.savefig("0Vs1_specplotcompare_col3.png")
plt.close()
"""
fig2=plt.figure()
plt.plot(f1, s1n, color="red", label="new fft function col1")
plt.plot(f2, s2n, color="blue", label="old fft function col1")
#plt.ylim(-0.5e-18,0.5e-18)
plt.yscale('log')
plt.xscale('log')
plt.legend()
fig2.savefig("specplotcompare_col1.png")
plt.close()

fig3=plt.figure()
plt.plot(f1, s1m, color="red", label="new fft function col2")
plt.plot(f2, s2m, color="blue", label="old fft function col2")
#plt.ylim(-0.5e-18,0.5e-18)
plt.yscale('log')
plt.xscale('log')
plt.legend()
fig3.savefig("specplotcompare_col2.png")
"""

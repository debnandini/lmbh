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

"""def plots(file):
	file = str(file)
	ft, A, E, T = load_text(file)
	c = "dgb"
	#ftc, Ac, Ec, Tc = load_text("file.split('.')[0]+"_dgb.dat")
	ftc, Ac, Ec, Tc = load_text("./seg_%s/"%(c)+file)
	At = A-Ac
	Et = E-Ec
	Tt = T-Tc
	 ftt = np.delete(ft, np.nonzero(np.isin(ft, ftc, assume_unique=False, invert=False)))
	At = np.delete(A, np.nonzero(np.isin(A, Ac, assume_unique=False, invert=False)))
	Et = np.delete(E, np.nonzero(np.isin(E, Ec, assume_unique=False, invert=False)))
	Tt = np.delete(T, np.nonzero(np.isin(T, Tc, assume_unique=False, invert=False)))
	print("len(At), len(ftt), len(Et), len(Tt)", len(At), len(ftt), len(Et), len(Tt))
	print("len(A), len(ft), len(E), len(T)", len(A), len(ft), len(E), len(T))
	print("len(Ac), len(ftc), len(Ec), len(Tc)", len(Ac), len(ftc), len(Ec), len(Tc))
	print("len(np.nonzero(np.isin(ft, ftc, assume_unique=False, invert=False)))", len(np.nonzero(np.isin(ft, ftc, assume_unique=False, invert=False))[0]), np.nonzero(np.isin(ft, ftc, assume_unique=False, invert=False))[0])
	print("len(np.nonzero(np.isin(A, Ac, assume_unique=False, invert=False)))", len(np.nonzero(np.isin(A, Ac, assume_unique=False, invert=False))[0]), np.nonzero(np.isin(A, Ac, assume_unique=False, invert=False))[0])
	print("len(np.nonzero(np.isin(E, Ec, assume_unique=False, invert=False)))", len(np.nonzero(np.isin(E, Ec, assume_unique=False, invert=False))[0]), np.nonzero(np.isin(E, Ec, assume_unique=False, invert=False))[0])
	print("len(np.nonzero(np.isin(T, Tc, assume_unique=False, invert=False)))", len(np.nonzero(np.isin(T, Tc, assume_unique=False, invert=False))[0]), np.nonzero(np.isin(T, Tc, assume_unique=False, invert=False))[0])
	#print("len(np.setdiff1d(Ac, A, assume_unique=False))", len(np.setdiff1d(Ac, A, assume_unique=False)))
	np.savetxt(file.split('.')[0]+"cut_%s.dat"%(c), np.transpose([ft,At,Et,Tt]))
	#np.savetxt("cut.dat", np.transpose([ft,At,Et,Tt]))
if __name__ == '__main__':
	plots(glob.glob('./AET_seg*_f.dat'))
"""
def cut(c, file):
	ft, A, E, T = load_text(file)
	ftc, Ac, Ec, Tc = load_text("./seg_%s/"%(c)+file)
	At = A-Ac
	Et = E-Ec
	Tt = T-Tc
	return(ft, At, Et, Tt)

files = sys.argv
files = files[1:]
#tocut = np.array(["dgb","vgb","igb"])
tocut = np.array(["dgb","vgb","igb"])
for file in files:
	for c in tocut:
		ft, At, Et, Tt = cut(c, file)
		np.savetxt(file, np.transpose([ft,At,Et,Tt]), fmt='%.15e')
		print(c, "has been cut from", file)


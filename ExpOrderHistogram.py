import re
import os
import math
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

def CharacteristicFunction(k, cumulants):
	lnG = 0.
	for i in range(1, len(cumulants)+1):
		lnG += np.power(1.j * k, i) / math.factorial(i) * cumulants[i-1]
	return np.exp(lnG)

def Normalize(values):
	for i in range(1, len(values)):
		total = 0.0
		for j in range(len(values[i])):
			total += values[i][j]
		if total == 0.0:
			total = 1.0
		for j in range(len(values[i])):
			values[i][j] /= total
	return values

def ParseQuantities(filename, cols):
	with open(filename) as inputFile:
		data = inputFile.read()
		values = []
		for i in range(cols):
			values.append([])
		for line in data.split('\n'):
			l = line.split()
			for i in range(len(l)):
				values[i].append(float(l[i]))
	return values

suffix = '.single.task0001/run1.exporderhist.txt'
filelist = ['job-T0.085-V1.3', 'job-T0.06-V0.7', 'job-T0.5-V1.5']
exporderCumulants = {}
exporderCumulants['job-T0.085-V1.3'] = [23.406550, 37.08387190]
exporderCumulants['job-T0.06-V0.7'] = [13.06175194, 18.188561887906218]
exporderCumulants['job-T0.09-V1.6'] = [30.92744575, 49.44462259]
exporderCumulants['job-T0.09-V1.5'] = [27.85101426, 44.54162779]
exporderCumulants['job-T0.5-V1.5'] = [4.84830623, 8.31184134]

for i in range(len(filelist)):
	if (not os.path.exists(filelist[i]+suffix)):
		continue
	values = ParseQuantities(filelist[i]+suffix, 3)
	values = Normalize(values)
	Z = sum(values[1])
	print filelist[i] + " " + str(values[1][0] / Z)
	plt.figure(0)
	plt.subplot(2, 2, i)
	plt.title(filelist[i])
	plt.xlabel("k")
	plt.ylabel("#")
	z, = plt.plot(np.array(values[0]), np.array(values[1]), "b-", linewidth=2.0, label=r'$Z_{k}$')
	mean = exporderCumulants[filelist[i]][0]
	sigma = np.sqrt(exporderCumulants[filelist[i]][1])
	gauss = mlab.normpdf(np.array(values[0]), mean, sigma)
	ed, = plt.plot(np.array(values[0]), gauss, "r-", linewidth=2.0, label=r'$Z^{ED}_{k}$')
	G_k = [[], []]
	kmin = values[0][0]
	kmax = values[0][-1]
	N = int(kmax - kmin)
	for k in range(N):
		G_k[0].append(kmin + (kmax - kmin) * k / N)
		G_k[1].append(CharacteristicFunction(2.*np.pi*k/N, exporderCumulants[filelist[i]]) / N)
	P_n = np.fft.fft(G_k[1]).real
	P_n -= min(P_n)
	P_n /= sum(P_n)
	ed2, = plt.plot(np.array(G_k[0]), P_n, "g-", linewidth=2.0, label=r'$Z^{ED}_{k}$')
		
	plt.legend()
plt.show()
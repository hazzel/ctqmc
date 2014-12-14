import re
import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from ParseDataOutput import *

#V = 1.5, T = 0.5
#ED = [0.25, -0.13467517, 0.09141432, -0.10371007]
#V = 1.5, T = 0.2
#ED = [0.25000000, -0.13996493, 0.09537386, -0.11612933]
#V = 1.5, T = 0.2
#ED = [0.25000000, -0.13996493, 0.09537386, -0.11612933]
#V = 1.6, T = 0.25
#ED = [0.25000000, -0.14676379, 0.10445700, -0.12268807]
#V = 1.5, T = 0.09
ED = [0.25000000, -0.13925507, 0.09439008, -0.11540501]
#V = 1.4, T = 0.07
#ED = [0.25000000, -0.13342669, 0.08676040, -0.11000114]

filelist = ['out/job-T0.09-V1.5-2.task0001.out']
for i in range(len(filelist)):
	if (not os.path.exists(filelist[i])):
		continue
	elist = ParseEvalables(filelist[i])
	y = ArrangePlot(elist, "Correlations\[\d+\]")
	x = np.arange(0., len(y[0]), 1.0)
	plt.figure(i)
	plt.subplot(1, 2, 1)
	plt.title(filelist[i])
	plt.plot(x, np.array(y[0]), "b-", linewidth=2.0, label=r'$C^{MC}_{R}$')
	plt.errorbar(x, np.array(y[0]), yerr=np.array(y[1]))
	plt.plot(x, np.array(ED), "r-", linewidth=2.0, label=r'$C^{ED}_{R}$')
	plt.xlabel("R")
	plt.ylabel("C(R)")
	plt.ylim([-0.3,0.3])
	
	z = []
	for j in range(len(ED)):
		z.append(y[0][j] - ED[j])
	plt.subplot(1, 2, 2)
	plt.plot(x, z, "ko", label=r'$C^{MC}_{R} - C^{ED}_{R}$')
	plt.errorbar(x, z, yerr=np.array(y[1]))
	plt.xlabel("R")
	plt.ylabel("C(R)")
	plt.xlim([x[0]-0.1,x[-1]+0.1])
	
	plt.legend()
	PrintData(x, y, "C(R)", filelist[i])
plt.show()
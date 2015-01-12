import re
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from ParseDataOutput import *

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

L = ["2", "3", "4", "5", "6", "7"]
for l in range(len(L)):
	filelist = glob.glob("out/job-mkl-L" + L[l] + "-V2.0-hex.task*.out")
	if len(filelist) == 0:
		continue
	filelist.sort()
	x = []
	y = []
	yerr = []
	for i in range(len(filelist)):
		if (not os.path.exists(filelist[i])):
			continue
		plist = ParseParameters(filelist[i])
		elist = ParseEvalables(filelist[i])
		x.append(float(plist["T"]))
		y.append( ArrangePlot(elist, "Binder")[0][0] )
		yerr.append( ArrangePlot(elist, "Binder")[1][0] )
		y = [i for j, i in sorted(zip(x, y))]
		x.sort()
	
	plt.figure(1)
	plt.title(filelist[0].split("/")[1].split(".task")[0])
	plt.xlabel("T")
	plt.ylabel("B")
	plt.plot(np.array(x), np.array(y), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
	plt.errorbar(np.array(x), np.array(y), yerr=np.array(yerr), color=color_cycle[l])
	plt.legend()
plt.show()

import re
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from ParseDataOutput import *

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
L = ["6"]
for l in range(len(L)):
	filelist = glob.glob("out/job-L" + L[l] + "-T0.25.task*.out")[1::3]
	filelist.sort()
	for i in range(len(filelist)):
		if (not os.path.exists(filelist[i])):
			continue
		x = []
		y = []
		plist = ParseParameters(filelist[i])
		elist = ParseEvalables(filelist[i])
		y = ArrangePlot(elist, "Correlations\[\d+\]")
		x = np.arange(0., len(y[0]), 1.0)
		
		plt.figure(1)
		plt.title(filelist[0].split("/")[1].split(".task")[0])
		plt.xlabel("R")
		plt.ylabel("C(R)")
		plt.plot(np.array(x), np.array(y[0]), "-", color=color_cycle[i%len(color_cycle)], linewidth=2.0, label=r'V='+plist["V"])
		plt.errorbar(np.array(x), np.array(y[0]), yerr=np.array(y[1]), color=color_cycle[i%len(color_cycle)])
		plt.legend()
plt.show()

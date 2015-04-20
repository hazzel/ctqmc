import re
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from ParseDataOutput import *

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

L = ["2", "3", "4", "5", "6", "9", "12", "15"]
for l in range(len(L)):
	filelist = []
	#filelist.append(glob.glob("plot_rhom_V2.0-T*/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_rhom_V2.0/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V1.5/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V1.625/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V1.75/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V1.875/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_hex_V2.0/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V2.0-T*/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V2.25/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V2.5/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V3.0/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_hex_V3.5/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_hex_V2.0-T*/*L" + L[l] + "*.out"))
	for f in range(len(filelist)):
		if len(filelist[f]) == 0:
			continue
		filelist[f].sort()
		x = []
		y = []
		yerr = []
		for i in range(len(filelist[f])):
			if (not os.path.exists(filelist[f][i])):
				continue
			plist = ParseParameters(filelist[f][i])
			elist = ParseEvalables(filelist[f][i])
			x.append(float(plist["T"]))
			y.append( ArrangePlot(elist, "Binder")[0][0] )
			yerr.append( ArrangePlot(elist, "Binder")[1][0] )
			y = [i for j, i in sorted(zip(x, y))]
			x.sort()
		
		plt.figure(f)
		plt.title(filelist[f][0])
		plt.xlabel("T")
		plt.ylabel("B")
		plt.plot(np.array(x), np.array(y), "o", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
		plt.errorbar(np.array(x), np.array(y), yerr=np.array(yerr), color=color_cycle[l])
		plt.legend(loc=2)
plt.show()

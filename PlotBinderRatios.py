import re
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from ParseDataOutput import *

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'b']

L = ["3", "4", "5", "6", "9", "12", "15"]
for l in range(len(L)):
	filelist = []
	#filelist.append(glob.glob("plot_rhom_V2.0-T*/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_rhom_V2.0/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_hex_V1.355/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_hex_V1.5/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_hex_V1.625/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_hex_V1.75/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_hex_V1.875/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_hex_V2.0/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_hex_V2.0-T*/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot/plot_hex_V2.25/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_hex_V2.5/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_hex_V3.0/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_hex_V3.5/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_hex_V4.0/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_hex_T0.08/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot/plot_rhom_V1.5/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot/plot_rhom_V1.625/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot/plot_rhom_V1.75/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_rhom_V1.8125/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_rhom_V1.9375/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot/plot_rhom_V1.875/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_rhom_V2.0/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot/plot_rhom_V2.125/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot/plot_rhom_V2.5/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot/plot_rhom_V3.0/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot/plot_rhom_V1.355/*L" + L[l] + "*.out"))
	flag = "m2"
	ze = 1.45
	
	for f in range(len(filelist)):
		if len(filelist[f]) == 0:
			continue
		filelist[f].sort()
		x = []
		yB = []
		yBerr = []
		yM2 = []
		yM2err = []
		ycM2 = []
		ycM2err = []
		for i in range(len(filelist[f])):
			if (not os.path.exists(filelist[f][i])):
				continue
			plist = ParseParameters(filelist[f][i])[0]
			elist = ParseEvalables(filelist[f][i])[0]
			#x.append((float(plist["T"]) - 0.46)/0.46 * float(L[l]))
			if flag == "c":
				x.append(float(plist["T"]) * float(L[l]))
				#if float(plist["T"]) > 0.25: continue
			else:
				x.append(float(plist["T"]))
			yB.append( ArrangePlot(elist, "Binder")[0][0] )
			yBerr.append( ArrangePlot(elist, "Binder")[1][0] )
			yM2.append( ArrangePlot(elist, "M2")[0][0] * float(L[l])**0.)
			yM2err.append( ArrangePlot(elist, "M2")[1][0] * float(L[l])**0.)
			ycM2.append( ArrangePlot(elist, "M2")[0][0] * float(L[l])**ze)
			ycM2err.append( ArrangePlot(elist, "M2")[1][0] * float(L[l])**ze)
			yB = [i for j, i, k in sorted(zip(x, yB, yBerr))]
			yBerr = [k for j, i, k in sorted(zip(x, yB, yBerr))]
			yM2 = [i for j, i, k in sorted(zip(x, yM2, yM2err))]
			yM2err = [k for j, i, k in sorted(zip(x, yM2, yM2err))]
			ycM2 = [i for j, i, k in sorted(zip(x, ycM2, ycM2err))]
			ycM2err = [k for j, i, k in sorted(zip(x, ycM2, ycM2err))]
			x.sort()
		
		#plt.xlim([0.25,0.36])
		plt.figure(f)
		plt.title(r"$t=1,\ V=" + plist["V"] + "$")
		#plt.title(filelist[f][0])
		#plt.subplot(1, 2, 1)
		if flag == "b":
			plt.xlabel(r"$T$")
			plt.ylabel(r"$B$")
			plt.plot(np.array(x), np.array(yB), "o", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l], markersize=5)
			plt.errorbar(np.array(x), np.array(yB), yerr=np.array(yBerr), color=color_cycle[l], markersize=5)
			plt.legend(loc=2)
		#plt.subplot(1, 2, 2)
		elif flag == "m2":
			plt.xlabel(r"$T$")
			plt.ylabel(r"$M_2 L^{\eta^{2D}}$")
			plt.plot(np.array(x), np.array(yM2), "o", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l], markersize=5)
			plt.errorbar(np.array(x), np.array(yM2), yerr=np.array(yM2err), color=color_cycle[l], markersize=5)
			plt.legend(loc=1)
		elif flag == "c":
			plt.title(r"$z+\eta=" + str(ze) + "$")
			plt.xlabel(r"$T L^z$")
			plt.ylabel(r"$M_2 L^{z+\eta}$")
			plt.plot(np.array(x), np.array(ycM2), "o", linestyle="None", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l], markersize=5)
			plt.errorbar(np.array(x), np.array(ycM2), yerr=np.array(yM2err), linestyle="None", color=color_cycle[l], markersize=5)
			plt.legend(loc=1)
plt.show()

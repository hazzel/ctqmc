import re
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from ParseDataOutput import *

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

L = ["2", "3", "4", "5", "6", "9", "12"]
for l in range(len(L)):
	filelist = []
	#filelist.append(glob.glob("plot_rhom_V2.0-T*/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_rhom_V2.0/*L" + L[l] + "*.out"))
	#filelist.append(glob.glob("plot_hex_V2.0/*L" + L[l] + "*.out"))
	filelist.append(glob.glob("plot_hex_V2.0-T*/*L" + L[l] + "*.out"))

	for f in range(len(filelist)):
		if len(filelist[f]) == 0:
			continue
		filelist[f].sort()
		z = 0
		beta = 0.125
		eta = 0.25
		nu = 1.
		V = 2.0
		Tc = 0.52
		x = []
		yM2 = []
		yM2err = []
		yM4 = []
		yM4err = []
		for i in range(len(filelist[f])):
			if (not os.path.exists(filelist[f][i])):
				continue
			plist = ParseParameters(filelist[f][i])
			elist = ParseEvalables(filelist[f][i])
			x.append(float(plist["T"]))
			yM2.append( ArrangePlot(elist, "M2")[0][0] * float(L[l])**(0.666) )
			yM2err.append( ArrangePlot(elist, "M2")[1][0] * float(L[l])**(0.666) )
			yM4.append( ArrangePlot(elist, "M4")[0][0] * float(L[l])**(2.*0.666) )
			yM4err.append( ArrangePlot(elist, "M4")[1][0] * float(L[l])**(2.*0.666) )
		
		plt.figure(f)
		m = re.search("V" + fpn, filelist[f][0])
		if m:
			plt.suptitle(r'$V = ' + m.group(0)[1:] + ',\ T_c = ' + str(Tc) + ',\ \\nu = ' + str(nu) + ',\ \\eta = ' + str(eta) + ',\ z = ' + str(z) + '$', fontsize=16)
		plt.subplot(2, 2, 1)
		plt.xlabel(r'$T$')
		plt.ylabel(r'$M_2 L^{z+\eta}$')
		plt.plot(np.array(x), np.array(yM2), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
		plt.errorbar(np.array(x), np.array(yM2), yerr=np.array(yM2err), color=color_cycle[l])
		plt.legend(loc='upper right')
		
		plt.subplot(2, 2, 3)
		plt.xlabel(r'$T$')
		plt.ylabel(r'$M_4 L^{2z+2\eta}$')
		plt.plot(np.array(x), np.array(yM4), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
		plt.errorbar(np.array(x), np.array(yM4), yerr=np.array(yM4err), color=color_cycle[l])
		plt.legend(loc='upper right')
		
		plt.subplot(2, 2, 2)
		plt.xlabel(r'$(T-T_c)L^{1/\nu}$')
		plt.ylabel(r'$M_2 L^{z+\eta}$')
		for i in range(len(x)):
			x[i] = (x[i] - Tc) * float(L[l])**(1./nu)
		plt.plot(np.array(x), np.array(yM2), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
		plt.errorbar(np.array(x), np.array(yM2), yerr=np.array(yM2err), color=color_cycle[l])
		plt.legend(loc='upper right')
		
		plt.subplot(2, 2, 4)
		plt.xlabel(r'$(T-T_c)L^{1/\nu}$')
		plt.ylabel(r'$M_4 L^{2z+2\eta}$')
		plt.plot(np.array(x), np.array(yM4), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
		plt.errorbar(np.array(x), np.array(yM4), yerr=np.array(yM4err), color=color_cycle[l])
		plt.legend(loc='upper right')
plt.show()

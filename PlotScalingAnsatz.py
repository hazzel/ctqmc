import re
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt
from ParseDataOutput import *

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k']

L = ["4", "5", "6", "7"]
for l in range(len(L)):
	filelist = glob.glob("out/job-L" + L[l] + "-T0.25.task*.out")
	if len(filelist) == 0:
		continue
	filelist.sort()
	z = 0
	eta = 0.25
	nu = 1.
	#T=0.3
	Vc = 1.85
	#T=0.25
	#Vc = 1.615
	x = []
	yM2 = []
	yM2err = []
	yM4 = []
	yM4err = []
	for i in range(len(filelist)):
		if (not os.path.exists(filelist[i])):
			continue
		plist = ParseParameters(filelist[i])
		elist = ParseEvalables(filelist[i])
		x.append(float(plist["V"]))
		yM2.append( ArrangePlot(elist, "M2")[0][0] * float(L[l])**(z+eta) )
		yM2err.append( ArrangePlot(elist, "M2")[1][0] * float(L[l])**(z+eta) )
		yM4.append( ArrangePlot(elist, "M4")[0][0] * float(L[l])**(2.0*(z+eta)) )
		yM4err.append( ArrangePlot(elist, "M4")[1][0] * float(L[l])**(2.0*(z+eta)) )
	
	plt.figure(1)
	m = re.search("T" + fpn, filelist[0])
	if m:
		plt.suptitle(r'$T = ' + m.group(0)[1:] + ',\ V_c = ' + str(Vc) + ',\ \\nu = ' + str(nu) + ',\ \\eta = ' + str(eta) + ',\ z = ' + str(z) + '$', fontsize=16)
	plt.subplot(2, 2, 1)
	plt.xlabel(r'$V$')
	plt.ylabel(r'$M_2 L^{z+\eta}$')
	plt.plot(np.array(x), np.array(yM2), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
	plt.errorbar(np.array(x), np.array(yM2), yerr=np.array(yM2err), color=color_cycle[l])
	plt.legend(loc='upper left')
	
	plt.subplot(2, 2, 3)
	plt.xlabel(r'$V$')
	plt.ylabel(r'$M_4 L^{2z+2\eta}$')
	plt.plot(np.array(x), np.array(yM4), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
	plt.errorbar(np.array(x), np.array(yM4), yerr=np.array(yM4err), color=color_cycle[l])
	plt.legend(loc='upper left')
	
	if l > 0:
		plt.subplot(2, 2, 2)
		plt.xlabel(r'$(V-V_c)L^{1/\nu}$')
		plt.ylabel(r'$M_2 L^{z+\eta}$')
		for i in range(len(x)):
			x[i] = (x[i] - Vc) * float(L[l])**(1./nu)
		plt.plot(np.array(x), np.array(yM2), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
		plt.errorbar(np.array(x), np.array(yM2), yerr=np.array(yM2err), color=color_cycle[l])
		plt.legend(loc='upper left')
		
		plt.subplot(2, 2, 4)
		plt.xlabel(r'$(V-V_c)L^{1/\nu}$')
		plt.ylabel(r'$M_4 L^{2z+2\eta}$')
		plt.plot(np.array(x), np.array(yM4), "-", color=color_cycle[l], linewidth=2.0, label=r'L='+L[l])
		plt.errorbar(np.array(x), np.array(yM4), yerr=np.array(yM4err), color=color_cycle[l])
		plt.legend(loc='upper left')
plt.show()

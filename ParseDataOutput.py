import re
import os
import numpy as np
import matplotlib
matplotlib.use('TkAgg') 
import matplotlib.pyplot as plt

fpn = '[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'
NameLine = 'Name\s+=\s+.+\r?\n'
BinsLine = 'Bins\s+=\s+.+\r?\n'
BinLengthLine = 'BinLength\s+=\s+.+\r?\n'
MeanLine = 'Mean\s+=\s+.+\r?\n'
ErrorLine = 'Error\s+=\s+.+\r?\n'
VarianceLine = 'Variance\s+=\s+.+\r?\n'
AutocorrelationtimeLine = 'Autocorrelationtime\s+=\s+.+\r?\n'

RebinLength = 'ReBinLength\s+=\s+' + fpn + '\s'
RebinBins = 'Bins\s+=\s' + fpn + '\s+'
RebinMean = 'Mean\s+=\s' + fpn + '\s+'
RebinError = 'Error\s+=\s' + fpn + '\s+'

Parameterset = 'Parameters'
QuantityBlock = NameLine + BinsLine + BinLengthLine + MeanLine + ErrorLine + VarianceLine + AutocorrelationtimeLine
EvalableBlock = NameLine + BinsLine + BinLengthLine + MeanLine + ErrorLine + VarianceLine
RebinLine = RebinLength# + RebinBins + RebinMean + RebinError

def ParseParameterBlocks(filename):
	with open(filename) as inputFile:
		data = inputFile.read()
		return data.split("PARAMETERS\n")[1:]

def ParseQuantityValue(block, lineexp):
	value = ''
	m = re.search(lineexp, block)
	if m:
		value = m.group(0).split('=')[1].strip()
	return value

def ParseParameters(filename):
	plist = []
	data = ParseParameterBlocks(filename)
	for block in data:
		p = {}
		for line in block.split("\n"):
			if "=" in line:
				name = line.split("=")[0].strip()
				if name == "Name":
					plist.append(p)
					break
				value = line.split("=")[1].strip()
				p[name] = value
	return plist
	
def ParseQuantities(filename):
	plist = []
	data = ParseParameterBlocks(filename)
	for block in data:
		subblockList = re.findall(QuantityBlock, block)
		quantities = []
		for i in range(len(subblockList)):
			q = Quantity()
			q.Parse(subblockList[i])
			quantities.append(q)
		qlist.append(quantities)
	return qlist

def ParseEvalables(filename):
	elist = []
	data = ParseParameterBlocks(filename)
	for block in data:
		subblockList = re.findall(EvalableBlock, block)
		evalables = []
		for i in range(len(subblockList)):
			q = Evalable()
			q.Parse(subblockList[i])
			evalables.append(q)
		elist.append(evalables)
	return elist
	
class BinningAnalysis:
	Rebinlength = []
	Bins = []
	Mean = []
	Error = []

class Quantity:
	Name = ''
	Bins = 1
	BinLength = 1
	Mean = 0.
	Error = 0.
	Variance = 0.
	Autocorrelationtime = 0.
	
	def Parse(self, block):
		self.Name = ParseQuantityValue(block, NameLine)
		self.Bins = ParseQuantityValue(block, BinsLine)
		self.BinLength = ParseQuantityValue(block, BinLengthLine)
		self.Mean = ParseQuantityValue(block, MeanLine)
		self.Error = ParseQuantityValue(block, ErrorLine)
		self.Variance = ParseQuantityValue(block, VarianceLine)
		self.Autocorrelationtime = ParseQuantityValue(block, AutocorrelationtimeLine)
		
class Evalable:
	Name = ''
	Bins = 1
	BinLength = 1
	Mean = 0.
	Error = 0.
	Variance = 0.
	
	def Parse(self, block):
		self.Name = ParseQuantityValue(block, NameLine)
		self.Bins = ParseQuantityValue(block, BinsLine)
		self.BinLength = ParseQuantityValue(block, BinLengthLine)
		self.Mean = ParseQuantityValue(block, MeanLine)
		self.Error = ParseQuantityValue(block, ErrorLine)
		self.Variance = ParseQuantityValue(block, VarianceLine)
		
def ArrangePlot(quantities, regex):
	y = [[], []]
	for q in quantities:
		m = re.match(regex, q.Name)
		if m:
			y[0].append(float(q.Mean))
			y[1].append(float(q.Error))
	return y
	
def PrintData(x, y, name, filename):
	print (name + " in " + filename + ":")
	for i in range(len(x)):
		print ("{:.1f}".format(x[i]) + " : " + "{:.5f}".format(y[0][i]) + " +- " + "{:.5f}".format(y[1][i]))
	print ("")

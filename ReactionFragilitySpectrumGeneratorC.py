#!/usr/bin/python

#print("# Reaction Fragility Spectrum Generator")
#print("# Wersja 0.6")
#print("# Data 2018.06.15")
#print("# Autor: Mateusz Jedrzejewski <mateusz.jedrzejewski@pwr.edu.pl>")

import math
import numpy

K = None
F = None
R = None
s = 1
ksi = None
dksi = 0.05368

#Natoms = self.getNumberOfAtoms()
Natoms = 4                    # <----- edit it
Nrows = Natoms * 3
Nfreq = Nrows - 6
#Nksi = 86
Nksi = 171                    # <----- edit it -1 
IRC = []

f = open("CHOF_IRC_PO_gas.log", 'r')      # <----- edit it

line = f.readline()
while line != '':
	# path direction
	if 'calculation of the REVERSE path.' in line:
		s = -1
	if 'NET REACTION COORDINATE UP TO THIS POINT' in line:
		line=line.split()
		ksi=s*float(line[8])
		IRC.append([ksi,R,F,K])
	if ' Point Number:   1          Path Number:   1' in line:
		IRC.append([0.0,R,F,K])
	# Read coordinates
	if 'Input orientation:' in line:
		R = numpy.zeros((Natoms,3), numpy.float64) 
		# Header row
		line = f.readline()
		line = f.readline()
		line = f.readline()
		line = f.readline()
		for j in range(Natoms):
			data = f.readline().split()
			R[j,0],R[j,1],R[j,2]=float(data[3]),float(data[4]),float(data[5])
		# Convert from Ang to Bohr
		R *= 1.88971616463207
	# Read force
	if 'Forces (Hartrees/Bohr)' in line:
		F = numpy.zeros((Natoms,3), numpy.float64) 	
		# Header row
		line = f.readline()
		line = f.readline()
		for j in range(Natoms):
			data = f.readline().split()
			F[j,0],F[j,1],F[j,2]=float(data[2]),float(data[3]),float(data[4])
	# Read force constant matrix
	if 'Force constants in Cartesian coordinates:' in line:
		K = numpy.zeros((Nrows,Nrows), numpy.float64)
		for i in range(int(math.ceil(Nrows / 5.0))):
			# Header row
			line = f.readline()
			# Matrix element rows
			for j in range(i*5, Nrows):
				data = f.readline().split()
				for k in range(len(data)-1):
					K[j,i*5+k] = float(data[k+1].replace('D', 'E'))
					K[i*5+k,j] = K[j,i*5+k]
		# Convert from atomic units (Hartree/Bohr_radius^2) to J/m^2
		#K *= 4.35974417e-18 / 5.291772108e-11**2
		#IRC.append([0.0,R,F,K])
	# Read frequencies
	if 'Frequencies --' in line:
		Q =  []
		m = 0
		NfreqRows = Nfreq / 3
		for i in range(NfreqRows):
			data = line.split()
			for k in range(2,len(data)):
				Q.append(float(data[k]))
			for j in range(Natoms+7):
				line = f.readline()
	line = f.readline()
# Close file when finished
f.close()

IRC.sort()

def write_k_AA_csv():
	text_header = "IRC;k_AA_sum"
	for i in range(Natoms):
			text_header = text_header + ";k_AA_" + str(i+1)
	text_csv = text_header + "\n"
	for k in range(Nksi):
		k_AA=0
		C = numpy.zeros((Natoms,Natoms), numpy.float64)
		for i in range(Nrows):
			k_AA+=IRC[k][3][i][i]
		for i in range(Natoms):
			for j in range(Natoms):
				C[i,j]=0
				for m in range(3):
					C[i,j]+=IRC[k][3][3*i+m][3*j+m]
		k_AA=0
		for i in range(Nrows):
			k_AA += IRC[k][3][i][i]
		line = str(IRC[k][0])+";"+str(k_AA)+";"
		first = True
		for i in range(Natoms):
			if not first:
				line=line+";"
			line=line+str(C[i,i])
			first= False
		text_csv = text_csv + str(line) + "\n"
	return text_csv

#def derivative(x1,x2,x3):
#	return float(0.5*(x1+x2+x3))

def derivative(y1,y2,delta):
	return float((y1-y2)/(delta))

def write_a_ksi_csv():
	text_header = "IRC;a_ksi_sum"
	for i in range(Natoms):
		text_header = text_header + ";a_ksi_" + str(i+1)
	text_csv = text_header + "\n"
	C_tmp = numpy.zeros((Natoms,Natoms), numpy.float64)
	C = numpy.zeros((Natoms,Natoms), numpy.float64)
	for k in range(Nksi):
		C_tmp = C
		C = numpy.zeros((Natoms,Natoms), numpy.float64)
		a = numpy.zeros((Natoms,Natoms), numpy.float64)
		a_ksi=[]
		for i in range(Natoms):
			for j in range(Natoms):
				C[i,j]=0
				for m in range(3):
					C[i,j]+=IRC[k][3][3*i+m][3*j+m]
		a_ksi_sum=0
		for i in range(Natoms):
			a_tmp=derivative(C[i,i],C_tmp[i,i],dksi)
			a_ksi_sum+=a_tmp
			a_ksi.append([IRC[k][0],a_tmp])
		if k != 0:
			line = str(IRC[k-1][0])+";"+str(a_ksi_sum)+";"
			first = True
			for i in range(Natoms):
				if not first:
					line=line+";"
				line=line+str(a_ksi[i][1])
				first= False
			text_csv = text_csv + str(line) + "\n"
	return text_csv

def write_freq_ksi_csv():
	text_header = "IRC"
	for i in range(Nfreq):
			text_header = text_header + ";freq_ksi_" + str(i+1)
	text_csv = text_header + "\n"
	for k in range(Nksi):
		line = IRC[k][0]
		for i in range(Nfreq):
			line = line + ";" + str(IRC[k][4][i])
		text_csv = text_csv + str(line) + "\n"
	return text_csv

o = open("wyniki_k_AA.csv","w")
o.write(write_k_AA_csv())
o.close()
o = open("wyniki_a.csv","w")
o.write(write_a_ksi_csv())
o.close()
#o = open("wyniki_freq.csv","w")
#o.write(write_freq_ksi_csv())
#o.close()

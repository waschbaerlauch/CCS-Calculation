# import tons of stuff, this will take some time in the first run of the script
import sys

from scipy import optimize
from scipy import asarray as ar,exp
from scipy import stats
import scipy

from PyQt4.QtCore import SIGNAL
from PyQt4.QtGui import QFileDialog
from PyQt4 import QtGui, QtCore, uic

import matplotlib.pyplot as plt
import operator 

import numpy as np
import os, fnmatch
import subprocess

import csv
from numpy import *
import math
from pylab import *



# seperate thread to open the calculated results after running the script
class OpenSummary(QtCore.QThread):	
	def __init__(self):
		QtCore.QThread.__init__(self)
	
	def run(self):
		subprocess.Popen(UpdateThread.OpenOutput, shell=True).wait()
	
# seperate thread to calculate CCS from given input data. Calculation will not interfere with mainwindow
class UpdateThread(QtCore.QThread):	
	def __init__(self):
		QtCore.QThread.__init__(self)
	
	def run(self):
		# current path of the script
		rootdir = os.path.dirname(os.path.realpath(__file__))
		# create all needed lists and dictionaries
		results = {}			
		output = []
		T_list = []
		P_list = []
		datax = []	
		datay = []
		
		# call options stands for the string, which stands behind the CDCReader.exe and 
		# describes the resolving power of the translation of the measured data.
		# 'add' is a variable, which corresponds to the field of 'peak sensitivity' in the main window
		# Limit corresponds to the given limit in the main window
		add = str(myGUI.peak_sens)
		Limit = float(myGUI.peak_limit)		
		call_options = '--im_bin '+add	
		
		
		# iterate through "path" and find all files with the wanted ending, e.g. '*.txt.'
		# return all files with path+filename
		def findFiles (path, filter):
			for root, dirs, files in os.walk(path):
				for file in fnmatch.filter(files, filter):
					yield os.path.join(root, file)
	
		# definition of the gaussian fits
		# the multiple gaussian fit just describes the sum of many gaussians
		#it returns the calculated y-values
		def gaussian(x, height, centre, width):
			return height * np.exp( - (x - centre)**2.0 / (2.0 * width**2.0) )
		
		def two_gaussian(x, h1, c1, w1, h2, c2, w2):
			return gaussian(x, h1, c1, w1)+gaussian(x, h2, c2, w2)

		def three_gaussian(x, h1, c1, w1, h2, c2, w2, h3, c3, w3):
			return two_gaussian(x, h1, c1, w1, h2, c2, w2)+gaussian(x, h3, c3, w3)

		# open file to read and split the lines to make it possible to iterate through lines
		# returns readable lines
		def GetLines(InputFileName):
			InputFile = open(InputFileName,'r') 
			Lines = InputFile.read().splitlines() 
			InputFile.close()
			return Lines

		# searching for key words in textFile = _extern file and calculating the drift voltage
		def GetVoltage(textFile):
			for line in GetLines(textFile):
				if line.find('ADC Pusher Frequency') >= 0:
					PusherInterval = float(line[25:])/1000
					#print "time:" + str(PusherInterval)				
				if line.find('Helium Cell DC') >= 0:
					Helium_Cell_DC = float(line[15:])
					#print "He Cell DC:" + str(Helium_Cell_DC)			
				if line.find('Helium Exit') >= 0:
					Helium_Exit = float(line[15:])
					#print "Helium Exit:" + str(Helium_Exit)	
				if line.find('IMSBias') >= 0:
					IMS_Bias = float(line[10:])
					#print "IMS Bias:" + str(IMS_Bias)	
				if line.find('Transfer DC Entrance') >= 0:
					Transfer_DC_Entrance = float(line[22:])
					#print "Transfer DC Entrance:" + str(Transfer_DC_Entrance)
			Voltage = float((Helium_Cell_DC + Helium_Exit + IMS_Bias - Transfer_DC_Entrance)*(1.0-2.0/170.0))
			return Voltage
		
		# searching for measuring time and date in text file = _header.txt and returns a timestamp
		def GetDate(textFile):
			for line in GetLines(textFile):
				if line.find('$$ Acquired Date:') >= 0:
					date = line[18:]
					#print date			
				if line.find('$$ Acquired Time:') >= 0:
					time = line[18:]
					#print time
			return str(date+' '+time)
			
		# searching for matching times between measuring and Logfile
		# returning the assigned Temperature and Pressure for each timestamp
		def GetTemp(timesearch):
			data = np.genfromtxt(str(myGUI.logfile), names=True, dtype=("|S11", "|S5", float, "|S2", float))
			search_date = timesearch[:-9]
			search_time = timesearch[12:-3]			
			match = np.where((data['Date']==search_date)&(data['Time']==search_time))
			try:
				Temp = data[match][0][4]
				Driftgas = data[match][0][3]	
				Pressure = data[match][0][2]	
			except IndexError:
				Temp = 0
				Driftgas = 0
				Pressure = 0
			return Temp, Pressure, Driftgas

		# searching for key words in textFile = _extern file and returning the Pusher Interval
		# the factor of 0.00025 is empirical found and describes the mismatch between displayed pusher interval in the files and the real pusher interval
		def GetPusherInterval(textFile):
			for line in GetLines(textFile):
				if line.find('ADC Pusher Frequency') >= 0:
					PusherInterval_old = float(line[25:])/1000					
					PusherInterval = PusherInterval_old + 0.00025					
					#print "time:" + str(PusherInterval)
			return PusherInterval
				
		# the input x,y-values are to put in as arrays and corresponds to the drift times and intensities of the wanted peak
		# script will try to first fit single gaussian and see the residual error in % 
		# if error is <25%, it will output the value for a single gaussian, otherwise it will continue will multiple gaussian fits		
		# the optimated mean value will be returned
		def GetDriftTime(x, y):											
			try:	
				# this will find the biggest peak in data and use this as starting guess for the gauss fits
				index_max_int = np.nonzero(y == max(y))[0][0]
				total = y.sum()
				# this describes the starting 100 percent to evaluate the residual error after each run
				one_gauss = 100
				two_gauss = 100
				three_gauss = 100
				# this loop will try out several starting points in order to find the optimal parameters
				for i in range(int(x[index_max_int]-20),int(x[index_max_int]+20),2):
					errfunc1 = lambda p, x, y: (gaussian(x, *p) - y)**2
					initial_guess1 = [max(y),i,1]	
					opt_fit1, success1 = optimize.leastsq(errfunc1, initial_guess1[:], args =(x, y),maxfev=3000)
					err1 = np.sqrt(errfunc1(opt_fit1, x, y)).sum()/total*100
					# simply compare the results of each rounds and see whats the best fit
					# each new best fit will replace the old value
					if err1 < one_gauss:
						one_gauss = err1
						best_fit = opt_fit1
					else:		
						continue
				# if the best fit of single gaussian fit is good enough, it will stop here, else multiple gaussians are fitted
				if one_gauss < 22:
					fit = 'single gaussian'
				else:
					for i in range(int(x[index_max_int]-50),int(x[index_max_int]+50),2):
						errfunc2 = lambda p, x, y: (two_gaussian(x, *p) - y)**2
						initial_guess2 = [max(y),x[index_max_int],1,max(y)/3,i,1]
						opt_fit2, success2 = optimize.leastsq(errfunc2, initial_guess2[:], args =(x, y),maxfev=3000)
						err2 = np.sqrt(errfunc2(opt_fit2, x, y)).sum()/total*100
						if err2 < two_gauss:
							two_gauss = err2
							best_fit = opt_fit2
						else:
							continue	
					
					if two_gauss < 22:
						fit = 'double gaussian'
					
					else:
						for i in range(int(x[index_max_int]-50),int(x[index_max_int]+50),2):
							errfunc3 = lambda p, x, y: (three_gaussian(x, *p) - y)**2
							initial_guess3 = [max(y),x[index_max_int],1,max(y)/3,i,1, max(y)/3,i+10,1]
							opt_fit3, success3 = optimize.leastsq(errfunc3, initial_guess3[:], args =(x, y),maxfev=3000)
							err3 = np.sqrt(errfunc3(opt_fit3, x, y)).sum()/total*100
							if err3 < three_gauss:
								three_gauss = err3
								best_fit = opt_fit3								
							else:
								continue
						fit = 'triple gaussian'									
						
				# this is for the output of three drift times, even if there is only one available drift time
				if fit == 'single gaussian':
					time1 = best_fit[1]*PushInt
					time2 = False
					time3 = False	
					error = one_gauss
				elif fit == 'double gaussian':
					time1 = best_fit[1]*PushInt
					time2 = best_fit[4]*PushInt
					time3 = False
					error = two_gauss
				else:	
					time1 = best_fit[1]*PushInt
					time2 = best_fit[4]*PushInt
					time3 = best_fit[7]*PushInt
					error = three_gauss
				return time1, time2, time3, best_fit, error
			except ValueError:
				return False, False, False

		# voltage list x gets reversed for plotting, y-values are the drift times of the peaks
		# linregress plots a linear fit
		# Charge and Drift_Mass will be put in by the search list, T and P will be searched inside of the Logfile (GetTemp)	
		# sees if the results are not real numbers (nan) and then returns the calculated CCS and the error of this value
		def GetCCS(x, y):
			x2 = [1/value for value in x]	
			xaxis = ar(x2)
			yaxis = ar(y)	
			slope, intercept, r_value, p_value, std_err = scipy.stats.linregress(xaxis, yaxis)		
			K = (25.05**2)/(slope/1000)			
			K0 = K*(273.15/T)*(P/760)			
			Mass = content * Charge			
			reducedMass = (Mass*Drift_Mass)/(Mass+Drift_Mass)
			CCS = (Charge/K0)*((1/(reducedMass*T))**0.5)*18495.88486	
			error = (std_err/slope)*CCS
			error_CCS = (error/CCS)*100	
			if math.isnan(CCS) == True or math.isnan(error_CCS) == True:
				return str(False), str(False)
			else:
				return CCS, error_CCS

		# Truncates a float f to n decimal places without rounding
		def truncate(f, n):			
			s = '{}'.format(f)
			if 'e' in s or 'E' in s:
				return '{0:.{1}f}'.format(f, n)
			i, p, d = s.partition('.')
			return float('.'.join([i, (d+'0'*n)[:n]]))	
		
				
		# completed is the value for the progress bar in the mainwindow
		# the '1' is just for the user to see, that the calculation has started
		myGUI.completed = 1						
		self.emit(QtCore.SIGNAL('calc_progress'), myGUI.completed)				
		
		# list_of_raw is the list of data, which is selected by the user to get calculated			
		for raw in myGUI.list_of_raw:
			head, tail = os.path.split(raw)	
			#create output folder
			outputfolder_dir = 'Results_'+tail[:-6]
			output_dir = os.path.join(head, outputfolder_dir)				
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
									
			# find CDCReader.exe and perform conversion of data into text files
			for Exe in findFiles(rootdir, '*.exe'):
				Exe_path = os.path.join(rootdir, Exe)				
				subprocess.Popen(Exe_path+' '+raw+' '+call_options, shell=True).wait()
							
			# find _extern.inf file and get voltage and pusher interval out of it
			for extern in findFiles(raw, '*_extern.inf'):
				Voltage = GetVoltage(extern)	
				PushInt = GetPusherInterval(extern)		
				
			# find timestamp of the measured data in the _header.txt file
			for synaptLog in findFiles(raw, '*_HEADER.TXT'):
				timestamp = GetDate(synaptLog)	
				
			# calc will update the 'completed' value for the progress bar to inform the user, that one file is finished processing
			calc = 100/len(myGUI.list_of_raw)		
			myGUI.completed += calc						
			self.emit(QtCore.SIGNAL('calc_progress'), myGUI.completed)	
			
			# this will search for the temperature and pressure inside the logfile and append it to a temporary list
			if myGUI.parameter_checkbox.isChecked():
				temperature, pressure, Drift_Mass_new  = GetTemp(timestamp)
				if temperature == 0 or pressure == 0 or Drift_Mass_new == 0:
					continue
				else:
					T_list.append(temperature)
					P_list.append(pressure)

			# get data from current file
			data = np.genfromtxt(raw+"_1_imms.txt")	
			
			# remove useless ms.txt and imms.txt files which automatically get created while converting with CDCReader
			for too_much_files in findFiles(head, '*_1_ms.txt'):			
				os.remove(too_much_files)
			for too_much_files2 in findFiles(head, '*_1_imms.txt'):			
				os.remove(too_much_files2)
			

			# this gets the input searchlist as variable reader		
			reader = np.genfromtxt(str(myGUI.checkPeaklist))
			
			
			# see if the searchlist only has one peak (m/z and charge) in it or more
			if reader.size == 2:
				wanted_peak = reader[0]
				charge = reader[1]
				# then search for this peak in the measured data
				low_limit = wanted_peak-Limit
				high_limit = wanted_peak+Limit
				rows= np.where((data[:,0]>low_limit)&(data[:,0]<high_limit))
				# found peak is for the file output later
				found_peak = data[rows][0][0]
				# this is the found data to calculate CCS
				x = data[rows][:,1]
				y = data[rows][:,2]			
				Drift_Time1, Drift_Time2, Drift_Time3, best_fit,error = GetDriftTime(x, y)			
				# this multiple if statements look for single or multiple CCS (depend on the best fit of gaussian)
				if Drift_Time1 != False:
					if found_peak not in results:
						results[found_peak] = []
					results[found_peak].append([Voltage, Drift_Time1, charge])
					gauss="single"
					if Drift_Time2 != False:
						double_peak = found_peak+0.001
						if double_peak not in results:
							results[double_peak] = []
						results[double_peak].append([Voltage, Drift_Time2, charge])	
						gauss="double"
						if Drift_Time3 != False:
							triple_peak = found_peak+0.002
							if triple_peak not in results:
								results[triple_peak] = []
							# create the results dictionary to save the found data for later
							results[triple_peak].append([Voltage, Drift_Time3, charge])
							gauss="triple"
				#test image of gaussian fit
				if gauss == "single":
					self.emit(QtCore.SIGNAL('single'), best_fit,x,y, Voltage, wanted_peak, output_dir)
											
				if gauss == "double":	
					self.emit(QtCore.SIGNAL('double'), best_fit,x,y, Voltage, wanted_peak, output_dir)
					
				if gauss == "triple":
					self.emit(QtCore.SIGNAL('triple'), best_fit,x,y, Voltage, wanted_peak, output_dir)
			else:
				# this is doing the same like before, only this has more than one peak in the searchlist
				# it is not pretty, but it works :-)
				for wanted_peak, charge in reader:					
					low_limit = wanted_peak-Limit
					high_limit = wanted_peak+Limit
					rows= np.where((data[:,0]>low_limit)&(data[:,0]<high_limit))
					found_peak = data[rows][0][0]
					x = data[rows][:,1]
					y = data[rows][:,2]
					
					
					Drift_Time1, Drift_Time2, Drift_Time3, best_fit,error = GetDriftTime(x, y)			
					if Drift_Time1 != False:
						if found_peak not in results:
							results[found_peak] = []
						results[found_peak].append([Voltage, Drift_Time1, charge])
						gauss="single"
						
						if Drift_Time2 != False:
							double_peak = found_peak+0.0001
							if double_peak not in results:
								results[double_peak] = []
							results[double_peak].append([Voltage, Drift_Time2, charge])
							gauss="double"
							
							if Drift_Time3 != False:
								triple_peak = found_peak+0.0002
								if triple_peak not in results:
									results[triple_peak] = []
								results[triple_peak].append([Voltage, Drift_Time3, charge])	
								gauss="triple"
								
					
					#test image of gaussian fit
					if gauss == "single":
						self.emit(QtCore.SIGNAL('single'), best_fit,x,y, Voltage, wanted_peak, output_dir)
											
					if gauss == "double":	
						self.emit(QtCore.SIGNAL('double'), best_fit,x,y, Voltage, wanted_peak, output_dir)
					
					if gauss == "triple":
						self.emit(QtCore.SIGNAL('triple'), best_fit,x,y, Voltage, wanted_peak, output_dir)
					


				
		head, tail = os.path.split(raw)							
		
		# if progress bar has not reached 100 percent (strange update system), this will get it to 100 percent
		myGUI.completed = 100						
		self.emit(QtCore.SIGNAL('calc_progress'), myGUI.completed)
		
		if myGUI.parameter_checkbox.isChecked():			
			if Drift_Mass_new == 'N2':
				Drift_Mass = 28
			elif Drift_Mass_new == 'He':
				Drift_Mass = 4
			# get an average temperature and pressure out of the temporary lists
			T_arr = np.array(T_list)
			P_arr = np.array(P_list)
			T_new = np.mean(T_arr)
			if np.std(T_arr) > 1:
				for item in T_list:
					if abs(item-T_new)>3:							
						self.emit(QtCore.SIGNAL('warning'))
			P = np.mean(P_arr)
			T = 273.15 + T_new						
		else:			
			# check for the drift gas in the mainwindow
			if myGUI.nitrogen_button.isChecked():			
				Drift_Mass_new = 'N2'
				Drift_Mass = 28
			elif myGUI.helium_button.isChecked():			
				Drift_Mass_new = 'He'
				Drift_Mass = 4			
			T_new = float(myGUI.lineEdit.text())
			T = 273.15 + T_new	
			P = float(myGUI.lineEdit_2.text())		
		
		# iterate through results and append the x,y-values for single peaks in temporary list
		# calculate CCS out of temporary x,y-list and append to output list 
		# empty datax and datay list for next peak
		
		
		#Name = raw[:-6]+".2summary.csv"
		#OutPut = open(Name, 'w')
		
		#OutPut.write(results)
		#for item in results:		
			#for voltage, driftzeit, charge in results[item]:					
				#print item, voltage, driftzeit
		
		for content in results:		
			for vol, drift, Charge in results[content]:	
				#print vol,drift,Charge
				datax.append(vol)
				datay.append(drift)										
			CCS, error_CCS = GetCCS(datax, datay)
			if CCS != 'False' and error_CCS != 'False':
				output.append([content, CCS, error_CCS])
			datax = []	
			datay = []	
			
			
		# Sort list of results for m/z
		output_sorted = sorted(output, key=lambda x: (x[0]))	
		
		# sort searchlist for m/z
		if reader.size == 2:			
			mass_search = reader[0]
			charge_search = reader[1]
		else:
			search_output_sorted = sorted(reader, key=lambda x: (x[0]))
		
					
		# Write results and parameters to file.
		OutPutFileName = output_dir+'/'+tail[:-6]+".summary.csv"	
		
		# create name for open summary thread to call
		UpdateThread.OpenOutput = OutPutFileName
		# open new file to write (will always overwrite old files with same name)
		OutPutFile = open(OutPutFileName, 'w')
		# create csv format output (again two times, one for single peak search and one for multiple peak search)
		OutPutFile.write('Searched:,')
		OutPutFile.write(',')
		OutPutFile.write('m/z, Charge')
		OutPutFile.write('\n')
		if reader.size == 2:
			OutPutFile.write(',')
			OutPutFile.write(',')
			OutPutFile.write(str(mass_search,)+',')
			OutPutFile.write(str(charge_search,)+',')
			OutPutFile.write('\n')						
		else:
			for mass_search, charge_search in search_output_sorted:				
				OutPutFile.write(',')
				OutPutFile.write(',')
				OutPutFile.write(str(mass_search,)+',')
				OutPutFile.write(str(charge_search,)+',')
				OutPutFile.write('\n')

		OutPutFile.write('\n')
		OutPutFile.write('Found:,')
		OutPutFile.write(',')
		OutPutFile.write('m/z, CCS, Error in %,')
		OutPutFile.write('\n')
		for mass, CCS, error_CCS in output_sorted:
			mass_new = truncate(mass,2)
			CCS_new = format(CCS, '.2f')
			error_CCS_new = format(error_CCS, '.2f')
			OutPutFile.write(',')
			OutPutFile.write(',')
			OutPutFile.write(str(mass_new,)+',')
			OutPutFile.write(str(CCS_new,)+',')
			OutPutFile.write(str(error_CCS_new,)+',')
			OutPutFile.write('\n')
		
		T_new2 = format(T_new, '.2f')
		P_new = format(P, '.2f')
		OutPutFile.write('\n')
		OutPutFile.write('\n')
		OutPutFile.write('Parameter:,')	
		OutPutFile.write(',')
		OutPutFile.write('Drift Gas, P (Torr), T (C),')	
		OutPutFile.write('\n')
		OutPutFile.write(',')
		OutPutFile.write(',')
		OutPutFile.write(str(Drift_Mass_new,)+',')
		OutPutFile.write(str(P_new,)+',')
		OutPutFile.write(str(T_new,)+',')

		# Close the output file.
		OutPutFile.flush()
		OutPutFile.close()				
		
# this is the maindow and the graphical user interface
class myGUI(QtGui.QMainWindow):
	def __init__(self, parent=None):
		# start GUI
		super(myGUI, self).__init__()
		# load design
		uic.loadUi('myUI_new.ui', self)		
		# connect buttons to a function call
		self.logfile_check.setDisabled(True)
		self.select_logfile.setDisabled(True)
		self.connect(self.data_button,SIGNAL("clicked()"), self.select_data)
		self.connect(self.peaklist_button,SIGNAL("clicked()"), self.select_peaklist)		
		self.connect(self.commandLinkButton,SIGNAL("clicked()"), self.Check_before_run)
		self.connect(self.lineEdit, SIGNAL("editingFinished()"), self.Temperature)
		self.connect(self.lineEdit_2, SIGNAL("editingFinished()"), self.Pressure)
		self.connect(self.clear_input,SIGNAL("clicked()"), self.Clear_All)
		self.connect(self.open_summary,SIGNAL("clicked()"), self.Summary)
		self.connect(self.exitbtn,SIGNAL("clicked()"), self.Closing_Command)			
		self.connect(self.default_btn,SIGNAL("clicked()"), self.setBack)
		self.connect(self.change_btn,SIGNAL("clicked()"), self.setNew)
		self.connect(self.select_logfile,SIGNAL("clicked()"), self.Logfile)	
		self.connect(self.parameter_checkbox,SIGNAL("clicked()"), self.Change_look)	
		self.x = 0
		myGUI.peak_sens = self.peak_sens.text()	
		myGUI.peak_limit = self.peak_limit.text()
		myGUI.parameter_checkbox = self.parameter_checkbox
		
	# function calls, which the buttons are connected to
	# select to input temperature and pressure by hand or automatically
	def Change_look(self):	
			self.x += 1			
			if self.x%2==0:
				self.logfile_check.setAutoExclusive(False)
				self.logfile_check.setChecked(False)
				self.select_logfile.setDisabled(True)
			else:
				self.logfile_check.setAutoExclusive(False)
				self.logfile_check.setChecked(False)
				self.select_logfile.setDisabled(False)
	
	def gaussian(self, x, height, centre, width):
		return height * np.exp( - (x - centre)**2.0 / (2.0 * width**2.0) )
			
	#visualize
	def single_gauss(self, best_fit, x, y, voltage, peak, head):
		gauss1 = self.gaussian( x, *best_fit )
		plot( x, y, c='b', lw=3, label= 'measurement' )
		plot( x, gauss1, c='g',lw= 3, label='gaussian fit' )
		plt.axis([0, 200, 0, max(y)])
		plt.legend()
		output_dir	= head		
		sample_file_name = 'mass='+str(peak)+'_Vol='+str(int(voltage))+'_Fit=singlePeak.png'
		results_dir = os.path.join(output_dir, sample_file_name)		
		plt.savefig(results_dir)
		plt.clf()
		
		#file = open(name+'.txt', 'w+')		
		#xarray = np.array(x)
		#yarray = np.array(y)
		#table = np.array([xarray, yarray])
		#table = table.T
		#np.savetxt(file, table, fmt=['%d','%d'])
		#file.close()
		
	def double_gauss(self, data, x, y, voltage, peak, head):
		gauss1 = self.gaussian( x, data[0], data[1], data[2] )
		gauss2 = self.gaussian( x, data[3], data[4], data[5] )
		plot( x, y, c='b', lw=4, label= 'measurement' )
		plot( x, gauss1, c='g',lw= 3, label='first fit' )
		plot( x, gauss2, c='r',lw= 3, label='second fit' )
		plt.axis([0, 200, 0, max(y)])
		plt.legend()		
		output_dir	= head		
		sample_file_name = 'mass='+str(peak)+'_Vol='+str(int(voltage))+'_Fit=doublePeak.png'
		results_dir = os.path.join(output_dir, sample_file_name)		
		plt.savefig(results_dir)
		plt.clf()
		
		#file = open(name+'.txt', 'w+')		
		#xarray = np.array(x)
		#yarray = np.array(y)
		#table = np.array([xarray, yarray])
		#table = table.T
		#np.savetxt(file, table, fmt=['%d','%d'])
		#file.close()
	
		
	def triple_gauss(self, data, x, y, voltage, peak, head):
		gauss1 = self.gaussian( x, data[0], data[1], data[2] )
		gauss2 = self.gaussian( x, data[3], data[4], data[5] )
		gauss3 = self.gaussian( x, data[6], data[7], data[8] )
		plot( x, y, c='b', lw=4, label= 'measurement' )
		plot( x, gauss1, c='g',lw= 3, label='first fit' )
		plot( x, gauss2, c='r',lw= 3, label='second fit' )
		plot( x, gauss3, c='y',lw= 3, label='third fit' )
		plt.axis([0, 200, 0, max(y)])
		plt.legend()
		output_dir	= head		
		sample_file_name = 'mass='+str(peak)+'_Vol='+str(int(voltage))+'_Fit=triplePeak.png'
		results_dir = os.path.join(output_dir, sample_file_name)		
		plt.savefig(results_dir)
		plt.clf()
		
		#file = open(name+'.txt', 'w+')		
		#xarray = np.array(x)
		#yarray = np.array(y)
		#table = np.array([xarray, yarray])
		#table = table.T
		#np.savetxt(file, table, fmt=['%d','%d'])
		#file.close()
		
			
	# open file dialog to select logfile.txt
	def Logfile(self):
		myGUI.logfile = QtGui.QFileDialog.getOpenFileName(self,"Select Logfile", filter="*.txt")
		if myGUI.logfile != '':
			self.logfile_check.setEnabled(True)
			self.logfile_check.setChecked(True)		

	# this is for new entered peak sensitivity and peaks search limit
	def setNew(self):
		myGUI.peak_sens = self.peak_sens.text()	
		myGUI.peak_limit = self.peak_limit.text()	

	# put search sensitivity and limit back to default
	def setBack(self):
		self.peak_sens.setText('0.3')
		self.peak_limit.setText('0.15')	
		myGUI.peak_sens = self.peak_sens.text()	
		myGUI.peak_limit = self.peak_limit.text()
	
	# this updates the progress bar
	def progress(self, int):
		self.progressBar.setValue(int)
	
	# this shows a warning of too big temperature differences
	def update_warning(self):
		self.errorlabel.setText('Warning!')
		self.errorlabel_2.setText('One temperature variates more than')
		self.errorlabel_3.setText('three degrees celsius from average value!')
		
	# close mainwindow
	def Closing_Command(self):
		self.close()
	
	# open the summary of the last calculated results (available after first run)
	def Summary(self):
		self.OpenSummary = OpenSummary()		
		self.OpenSummary.start()		
		
	# clear all fields (like new run)
	def Clear_All(self):
		self.lineEdit.setText('')
		self.lineEdit_2.setText('')
		self.lineEdit_3.setText('')
		self.lineEdit_4.setText('')
		self.progressBar.setValue(0)
		self.logfile_check.setAutoExclusive(False)
		self.logfile_check.setChecked(False)
		self.logfile_check.setDisabled(True)
		self.errorlabel.setText('')
		self.errorlabel_2.setText('')
		self.errorlabel_3.setText('')
		
	# open file dialog to select multiple directories (the data set)
	def select_data(self):
		dialog = FileDialog(self)
		try:
			self.lineEdit_4.setText(str(len(dialog.raw_data)))
			myGUI.list_of_raw = dialog.raw_data			
		except AttributeError:
			self.lineEdit_4.setText('False')
		
	# open file dialog to select the search list with the wanted peaks in it (m/z and charge)
	def select_peaklist(self):
		adress = QtGui.QFileDialog.getOpenFileName(self,"Select Peaklist",filter="*.txt")		
		head, tail = os.path.split(str(adress))
		if adress == '':
			self.lineEdit_3.setText('Please enter valid peaklist')	
		else:
			self.lineEdit_3.setText(tail)		
		myGUI.checkPeaklist = adress	
	
	# field for temperature
	def Temperature(self):
		text = self.lineEdit.text()		
	
	# field for pressure
	def Pressure(self):
		text2 = self.lineEdit_2.text()

	# if 'run' is pressed, the script will first check the entered values
	def Check_before_run(self):	
		check_data = self.lineEdit_4.text()
		if check_data == '':
			self.lineEdit_4.setText('False')
			d = False
		elif check_data == 'False':
			d = False
		else:
			d = True
		check_peak = self.lineEdit_3.text()	
		if check_peak == '':
			self.lineEdit_3.setText('Please enter valid peaklist')
			e = False
		elif check_peak == 'Please enter valid peaklist':
			e = False
		else:
			e = True
		if myGUI.parameter_checkbox.isChecked() == False:			
			try:
				temp = float(self.lineEdit.text())
				self.errorlabel.setText('')
				a = True
			except ValueError:
				self.errorlabel.setText('Check the temperature!')
				a = False			
			try:
				press = float(self.lineEdit_2.text())	
				self.errorlabel_2.setText('')
				b = True
			except ValueError:
				self.errorlabel_2.setText('Check the pressure!')
				b = False		
			if self.nitrogen_button.isChecked() == True or self.helium_button.isChecked() == True:		
				c = True
				self.errorlabel_3.setText('')
			else:
				self.errorlabel_3.setText('Check the drift gas!')
				c = False							
				
			# if all entered values seem to be ok, the script will continue running the main calculation
			if a == True and b == True and c == True and d == True and e == True:				
				self.completed = 0
				myGUI.completed = self.completed
				myGUI.lineEdit_3 =self.lineEdit_3
				myGUI.nitrogen_button = self.nitrogen_button
				myGUI.helium_button = self.helium_button
				myGUI.lineEdit = self.lineEdit
				myGUI.lineEdit_2 = self.lineEdit_2
				self.CCS_Calculation()	
		else:			
			if self.x%2!=0:	
				if self.logfile_check.isChecked() == False:						
					f = False
					self.errorlabel.setText('Please enter valid Logfile!')
				else:
					f = True
					self.errorlabel.setText('')
			if d == True and e == True and f == True:
				self.errorlabel.setText('')
				self.errorlabel_2.setText('')
				self.errorlabel_3.setText('')
				self.CCS_Calculation()	
	# this is the main calculation and will be run by a seperate thread to prevent freezing and stuff
	def CCS_Calculation(self):
		self.UpdateThread = UpdateThread()		
		self.connect(self.UpdateThread,SIGNAL("calc_progress"), self.progress)	
		self.connect(self.UpdateThread,SIGNAL("warning"), self.update_warning)	
		self.connect(self.UpdateThread,SIGNAL("single"), self.single_gauss)
		self.connect(self.UpdateThread,SIGNAL("double"), self.double_gauss)
		self.connect(self.UpdateThread,SIGNAL("triple"), self.triple_gauss)
		self.UpdateThread.start()	

# this is the file dialog to select multiple directories
# it looks complicated because QWidget FileDialog is used to choose multiple files, not multiple directories like we need
class FileDialog(QtGui.QFileDialog):
	def __init__(dialog, parent=myGUI):
		QtGui.QFileDialog.__init__(dialog)
		dialog.setFileMode(dialog.Directory)
		dialog.setOption(dialog.ShowDirsOnly, True)		
		dialog.setFilter("Raw Data (*.raw)")		
		for view in dialog.findChildren((QtGui.QListView, QtGui.QTreeView)):
			if isinstance(view.model(), QtGui.QFileSystemModel):
				view.setSelectionMode(QtGui.QAbstractItemView.MultiSelection)		
		#this will create the list of data files we will process later
		if dialog.exec_():			
			dialog.filenames = dialog.selectedFiles()			
			dialog.raw_data = map(str, dialog.filenames)
			# this solves a problem of other path entries in the list than the data directories .raw we need
			for fnames in dialog.raw_data:
				if fnames[-4:] != '.raw':
					dialog.raw_data.remove(fnames)					
		dialog.show()	
		
def main():
	app = QtGui.QApplication(sys.argv)
	myapp = myGUI()
	myapp.show()
	app.exec_()
	
if __name__ == "__main__":
	main()






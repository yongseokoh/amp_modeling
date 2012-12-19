#Python function for SYSTOR'12, Analytic Modeling of SSD Write Performance 
#Original matlab function has been written by Peter Desnoyers, Northeastern University, pjd@ccs.neu.edu
#This was ported by Yongseok Oh (ysoh@uos.ac.kr), University of Seoul

# Copyright 2012 Yongseok Oh
# This file is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# It is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details. 
# You should have received a copy of the GNU General Public License
# along with this file. If not, see http://www.gnu.org/licenses/. 

from scipy.special import lambertw 
from scipy.optimize import fsolve 
from scipy.optimize import minimize_scalar
from boomslang import *
import sys
import math
import numpy as np


# Simple function for translating sf to alpha
def alpha(sf):
	return 1.0 / (1.0 - sf)

# Simple function for translating alpha to sf
def sf(alpha):
	return (alpha - 1.0)/alpha

# Eq.4
def lru(alpha): 
	t = alpha + lambertw(-alpha * math.exp(-alpha))
	A = alpha / t 
	return A

# Eq. 12
def greedy(alpha, np): 
	k = 1.0 + 1.0 / (2.0 * np)
	A = lru(k * alpha) / k
	return A

# Eq. 14
def lru_hc( alpha, r, f ): 
	hot_term = lambda x: r / (math.exp((r/f)*(alpha/x))-1)
	cold_term = lambda x: (1-r) / (math.exp(((1-r)/(1-f))*(alpha/x))-1)
	func = lambda x: 1 + hot_term(x) + cold_term(x) - x
	A = fsolve ( func, 1)
	return A

# Eq. 16
def greedy_hc( alpha, np, r, f): 
	k = 1.0 + 1.0 / ( 2.0 * np)
	A = lru_hc( k * alpha, r, f) / k 
	return A

# Eq. 20
def greedy_hcp( p, alpha, np, r, f):  
	if p > 1 or  p < 0:
	#print "Invalid p value", p
		return 1000

	alpha_h = (p*(alpha-1) + f) / f                           # Eq. 21
	alpha_c = ((1-p)*(alpha-1) + (1 - f)) / (1 - f)           # Eq. 22
	A = r * greedy(alpha_h, np) + (1-r) * greedy(alpha_c, np) # Eq. 20
	return A

# Optimal Greedy Cleaning
def greedy_hc_opt( alpha, np, r, f):

	func = lambda p: greedy_hcp(p, alpha, np, r, f)
	A = minimize_scalar( func, method = 'brent')
#	print A.x, func(A.x)
	return func(A.x)

# Written by ysoh 
def greedy_hcp2( p, alpha, np, r, f):  
	if p > 1 or  p < 0:
	#print "Invalid p value", p
		return 1000

	f_hot = float(f)
	active_rate = 0.2
	f_cold_active = active_rate * (1 - f) 
	f_cold_inactive = (1 - active_rate) * (1 - f)

	print "f_hot", f_hot, "f_cold_active", f_cold_active
#	print f_cold_active

	alpha_h = (p*(alpha-1) + f_hot) / f_hot                           # Eq. 21
	alpha_c = ((1-p)*(alpha-1) + (1 - f_cold_active)) / (1 - f_cold_active)           # Eq. 22

	Ahot = r * greedy(alpha_h, np) # Eq. 20
	Acold = (1 -r) * greedy(alpha_c, np)
#print Ahot, Acold 
	return Ahot , Acold, Ahot + Acold

""""
# Eq. 20
def lru_hcp( p, alpha, r, f):  
	if p > 1 or  p < 0:
	#print "Invalid p value", p
		return 1000

	alpha_h = (p*(alpha-1) + f) / f                           # Eq. 21
	alpha_c = ((1-p)*(alpha-1) + (1 - f)) / (1 - f)           # Eq. 22
	A = r * lru(alpha_h) + (1-r) * lru(alpha_c) # Eq. 20
	return A

# Optimal LRU Cleaning 
def lru_hc_opt( alpha, r, f):

	func = lambda p: lru_hcp(p, alpha, r, f)
	A = minimize_scalar( func, method = 'brent')
#	print A.x, func(A.x)
	return func(A.x)
"""
# Written by Eunjae Lee 
def real_traffic( alpha, r, f, length):
	func_str=""
	R_str=""
	for i in range(length):
		temp="math.exp( (-1.0) * %s * %s / ( %s * x ) )" % (alpha,r[i],f[i])
		R_str = "%s * %s / ( 1.0 - %s)" % (temp,r[i],temp)
		func_str=func_str+"+"+R_str

	func_str="lambda x: "+func_str+" + 1.0 - x "
	func=eval(func_str)
	A = fsolve(func,1.0)
	return A

# boomslang plot line function
def draw_line(xval, yval, label, linestyle, marker, color, linewidth = 1):
	line = Line()
	line.xValues = xval
	line.yValues = yval
	line.label = label
	line.lineStyle =  linestyle
	line.color = color
	line.marker = marker
	line.lineWidth = linewidth

	return line

# Translation function: disk u to victim u 
# This function was introduced by Menon, J. A Performance Comparison of RAID-5 and Log-Structured Arrays.
# This function is used by Hylog(FAST2002), Janus-FTL(EMSOFT2010), and OP-FCL(FAST2012)
def udtou(ud):
	threshold = 1e-9
	#int i, max;

	#double u1, u2, u;
	#double disku1, disku2, disku;

	i = 0
	max_iter = 10000
	disku1 = u1 = float(0)
	disku2 = u2 = float(1)
	
	while u2 - u1 > threshold:
		u = (u1+u2)/2

		disku = (u - 1.0) / math.log(u)
		if ud > disku:
			u1 = u
			disku1 = disku;
		else:
			u2 = u
			disku2 = disku;
	   
		i += 1
		if i > max_iter:
			print "Cannot find solution: %f" %ud
			#exit(1);
		
	#//printf("%lf\n", (double)(u1+u2)/(double)2);
	return float(u1+u2)/float(2) 

# gererating zipfian distribution CDF function 
def get_zipf(alpha, sample = 100, window_width = 100, gen_num = 100, cdf = 1):

	zipf_dist = np.random.zipf(alpha, sample)
	zipf_dist %= window_width
	zipf_dist = np.sort(zipf_dist, axis = 0)
	f, b =  np.histogram(zipf_dist, bins = window_width)
	f = np.sort(f, axis = 0)
	f = f[::-1]
	zipf_dist = f

	zipf_dist.astype('float32')
	zipf_dist = zipf_dist[:gen_num]
	if cdf == 1:
		zipf_dist = zipf_dist.cumsum()/float(sample)
	else:
		zipf_dist = zipf_dist/float(sample)

	n = zipf_dist.tolist();
			
	return n
#generating uniform distribution CDF function
def get_uniform(sample = 100, window_width = 100, gen_num = 100):

	uniform_dist = np.random.random_integers(0, window_width, (sample, ))

	uniform_dist %= window_width
	uniform_dist = np.sort(uniform_dist, axis = 0)
	f, b =  np.histogram(uniform_dist, bins = window_width)
#print "debug",  f, b
	f = np.sort(f, axis = 0)
	f = f[::-1]
	uniform_dist = f

	uniform_dist.astype('float32')
	uniform_dist = uniform_dist[:gen_num]
	uniform_dist = uniform_dist.cumsum()/float(sample)

	n = uniform_dist.tolist();
	return n

def validation():

	print "Uniform LRU Cleaning Eq. 4"
	print "Sf 0.03", lru(alpha(0.03))
	print "Sf 0.07", lru(alpha(0.07))

	print ""
	print "Uniform Greedy Cleaning Eq. 12"
	print "Sf 0.03", greedy(alpha(0.03), 64)
	print "Sf 0.05", greedy(alpha(0.05), 64)
	print "Sf 0.07", greedy(alpha(0.07), 64)
	print "Sf 0.09", greedy(alpha(0.09), 64)
	print "Sf 0.11", greedy(alpha(0.11), 64)
	print "Sf 0.17", greedy(alpha(0.17), 64)


	print ""
	print "Non-Uniform LRU Cleaning Eq. 14"
	print lru_hc(alpha(0.03), 0.90, 0.05)
	print lru_hc(alpha(0.07), 0.80, 0.2)
	print lru_hc(alpha(0.07), 0.90, 0.05)
	print lru_hc(alpha(0.11), 0.80, 0.2)
	print lru_hc(alpha(0.11), 0.90, 0.05)
	print lru_hc(alpha(0.20), 0.80, 0.2)
	print lru_hc(alpha(0.20), 0.90, 0.05)

	print ""
	print "Non-Uniform Greedy Cleaning Eq. 16" 
	print greedy_hc(alpha(0.03), 32, 0.9, 0.05)
	print greedy_hc(alpha(0.07), 64, 0.9, 0.05)
	print greedy_hc(alpha(0.07), 128, 0.8, 0.20)
	print greedy_hc(alpha(0.11), 64, 0.9, 0.05)
	print greedy_hc(alpha(0.11), 32, 0.8, 0.20)
	print greedy_hc(alpha(0.20), 64, 0.9, 0.05)
	print greedy_hc(alpha(0.20), 128, 0.8, 0.20)

	print ""
	print "Optimal Greedy Hot Cold Eq. 20"
	print greedy_hc_opt(alpha(0.07), 64, 0.9, 0.05)
	print greedy_hc_opt(alpha(0.07), 128, 0.8, 0.20)
	print greedy_hc_opt(alpha(0.11), 32, 0.8, 0.20)
	print greedy_hc_opt(alpha(0.11), 64, 0.9, 0.05)
	print greedy_hc_opt(alpha(0.20), 64, 0.9, 0.05)
	print greedy_hc_opt(alpha(0.20), 128, 0.8, 0.20)

	""""
	print ""
	print "Optimal LRU Hot Cold"
	print lru_hc_opt(alpha(0.07), 0.9, 0.05)
	print lru_hc_opt(alpha(0.07), 0.8, 0.20)
	print lru_hc_opt(alpha(0.11), 0.8, 0.20)
	print lru_hc_opt(alpha(0.11), 0.9, 0.05)
	print lru_hc_opt(alpha(0.20), 0.9, 0.05)
	print lru_hc_opt(alpha(0.20), 0.8, 0.20)
	"""


def validation2():
	print ""
	print "greedy cleaning for uniform traffic"
	for i in range(10, 100, 10):
	#print float(i)/100
		print float(i)/100, greedy(alpha(1.0-float(i)/100), 64)

	print ""
	print "greedy cleaning for non-uniform traffic"
	for i in range(10, 100, 10):
		print float(i)/100, greedy_hc(alpha(1.0-float(i)/100), 32, 0.9, 0.3)

def plot_workload():

	xvalues = []
	for i in range(1, 101, 1):
		xvalues.append(float(i)/100)

	alpha = 1.1
	zipf_curve = []
	for i in range(1, 10, 2):
		zipf_curve.append(get_zipf(alpha, 10000, 100, 100))
		alpha += 0.2

	line1 = draw_line(xvalues, zipf_curve[0], "zipf(1.1)", "-", "",  "blue" );
#line2 = draw_line(xvalues, zipf_curve[1], "zipf(1.3)", "-", "",  "green" );
#	line3 = draw_line(xvalues, zipf_curve[2], "zipf(1.5)", "-", "",  "pink" );
#	line4 = draw_line(xvalues, zipf_curve[3], "zipf(1.7)", "-", "",  "yellow" );
#	line5 = draw_line(xvalues, zipf_curve[4], "zipf(1.9)", "-", "",  "black" );

	uniform_curve = get_uniform(10000, 100, 100)
	line6 = draw_line(xvalues, uniform_curve, "uniform", "-", "",  "red" );

	plot = Plot()
	plot.add(line1)
#	plot.add(line2)
#	plot.add(line3)
#	plot.add(line4)
#	plot.add(line5)
	plot.add(line6)

	plot.xLabel = "SLC Utilization"
	plot.yLabel = "Update Frequency(CDF)"

	plot.hasLegend()
	plot.legendLabelSize = 10
	plot.yLimits = (0.0, 1.0)
	plot.xLimits = (0.0, 1.0)
	plot.setDimensions(4, 3, 100)
	plot.save ("workload.png")

	return 

def plot_amplification():
	np = 256
	cap_spare = 2.0 	# GB
	cap_hot = 2.0 	# Physcal Hot Space e.g. SLC
	cap_cold = 16.0 # Physical Cold Space e.g. MLC

	ln_slc_wa = [] # ln model wa for hot 
	ln_tlc_wa = []
	ysoh_slc_wa = []
	ysoh_tlc_wa = []

	for i in range(10, 100, 20):
		sf_util = float(100-i)/100
		sf_hot = cap_spare * sf_util
		sf_cold = cap_spare * (1.0 - sf_util)
		
		u_hot = (cap_hot-sf_hot)/cap_hot
		u_cold = (cap_cold - sf_cold)/cap_cold

		# ln modeling 
		ln_slc_wa.append(1.0/(1.0-udtou(u_hot)))
		ln_tlc_wa.append(1.0/(1.0-udtou(u_cold)))

		# approximation real workload
		# SLC part 
		logical_data = cap_hot - sf_hot 
		inactive_data = logical_data * 0.1
		migration_data = logical_data * 0.1
		physical_data = cap_hot

		logical_data -= (inactive_data + migration_data )
		physical_data -= inactive_data

		print "L", logical_data, cap_hot - sf_hot, logical_data/(cap_hot-sf_hot)

		u_hot = (logical_data)/physical_data 
		ysoh_slc_wa.append(1.0/(1.0-udtou(u_hot)))

		# TLC part 
		logical_data = cap_cold - sf_cold
		inactive_data = logical_data * 0.3
		migration_data = logical_data * 0.15
		physical_data = cap_cold

		logical_data -= (inactive_data + migration_data )
		physical_data -= inactive_data

		u_cold = logical_data/physical_data 
		ysoh_tlc_wa.append(1.0/(1.0-udtou(u_cold)))



	if len(sys.argv) != 2:
		print " Invalid parameters ... \n"
		exit(0)

	print "Data File: " + sys.argv[1];
#print "Output File: " + sys.argv[2];

	lines = Utils.getLinesFromFile( sys.argv[1], "(\d+.\d+)	(\d+.\d+)	(\d+.\d+)	(\d+.\d+)	(\d+.\d+)	(\d+.\d+)")

	if len(lines) == 0 :
		print " Invalid Data Format " + sys.argv[1]
		exit(0)
	line1 = draw_line(lines[0].xValues, ln_slc_wa, "predicted(SLC)", "-", "",  "blue", 2 );
	line2 = draw_line(lines[1].xValues, lines[1].yValues, "measured(SLC)", "-", "",  "red", 2 );
	line3 = draw_line(lines[3].xValues, ysoh_slc_wa, "improved(SLC)", "-", "", "green", 2);

	line4 = draw_line(lines[0].xValues, ln_tlc_wa, "predicted(TLC)", "--", "",  "blue", 2 );
	line5 = draw_line(lines[1].xValues, lines[3].yValues, "measured(TLC)", "--", "",  "red", 2 );
	line6 = draw_line(lines[3].xValues, ysoh_tlc_wa, "improved(TLC)", "--", "", "green", 2);



	plot = Plot()
	plot.add(line1)
	plot.add(line2)
	plot.add(line3)
	
	plot.add(line4)
	plot.add(line5)
	plot.add(line6)


	plot.xLabel = "SLC Utilization"
	plot.yLabel = "Write Amplification"
#plot.xTickLables = [0.1, 0.3, 0.5, 0.7, 0.9]
	plot.hasLegend()
	plot.legendLabelSize = 10
	plot.yLimits = (0, 9)
	plot.xLimits = ( min(lines[0].xValues), max(lines[0].xValues))
	plot.xLimits = (0.0, 1.0)
	plot.setDimensions(4, 3, 300)
	plot.save ("amplification.pdf")

####################### Main Start ######################

#validation()
#validation2()

def plot_real():

	sf = [float(i)/100 for i in range(1, 30, 1)]
	f = [float(1)/100 for i in range(0, 100)]
	r = get_zipf(1.1, 10000, 100, 100, cdf = 0)

	real = []
	grdy = []
	grdy_hc = []
	ln = []

	for i in sf:
		real.append(real_traffic(alpha(i), r, f, len(r)))
		grdy.append( greedy(alpha(i), 64) )
		grdy_hc.append(greedy_hc(alpha(i), 64, 0.7, 0.3))
		ln.append(1.0/(1.0-udtou(1.0-i)))

	line1 = draw_line( sf, real, "real", "-", "",  "red" );
	line2 = draw_line( sf, grdy, "greedy", "-", "",  "blue" );
	line3 = draw_line( sf, grdy_hc, "greedy_hc", "-", "",  "green" );
	line4 = draw_line( sf, ln, "ln", "-", "",  "pink" );

	plot = Plot()
	plot.add(line1)
	plot.add(line2)
	plot.add(line3)
	plot.add(line4)

	plot.xLabel = "Spare Factor"
	plot.yLabel = "Write Amplification"
	plot.hasLegend()
	plot.legendLabelSize = 10
#	plot.yLimits = (0, 9)
#plot.xLimits = ( min(sf), max(sf))
	plot.setDimensions(4, 3, 100)
	plot.save ("real.pdf")


if __name__ == "__main__":

	validation()
	plot_real()
	plot_workload()
#	plot_amplification()
	print " EOP " 


#################### End Program ##########################

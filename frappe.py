#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Description: VCF summary
from optparse import OptionParser
import os
import shlex, subprocess

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2011, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
frappe_dir = '/lgc/programs/frappe1.1/'

parser = OptionParser()

parser.add_option('-p', help='ped', dest='pedfile') #, nargs=2
parser.add_option('-m', help='map', dest='mapfile') #, nargs=2
parser.add_option('-o', help='prefix', dest='prefix') #, nargs=2
(options, args) = parser.parse_args()
filename = ".".join(options.pedfile.split("/")[-1].split(".")[:-1])
print filename


command = '../plink-1.07-x86_64/plink --noweb --file %s --missing' % (filename)
os.system(command)
die()
#recode ped map to 12
command = '../plink-1.07-x86_64/plink --noweb --file %s --recode12 --output-missing-genotype 0 --geno 0.1 --out %s_12' % (filename, filename)
os.system(command)


#counting individuals and markers before running frappe
output = subprocess.check_output("wc -l %s_12.ped" % (filename), shell=True)
n_individuals = int(output.split()[0])
output = subprocess.check_output("wc -l %s_12.map" % (filename), shell=True)
n_markers = int(output.split()[0])
outfile = open("frappeconfig.txt", 'w')
outfile.writelines("""GenotypeFile  ="%s_12.ped" ## Mandatory genotype file name
MaxIter=  10000    ## maximum number of EM iterations
K      =  3  ## Number of ancestral individuals
M      = %s  ## Number of markers
I      =  %s  ## Number of individuals;
step   =  200 ### Define either step or Nout
Nout   =  0
printP  = 1     ## optional
threshold= 10000  ## optional convergence threshold

""" % (filename,n_markers,n_individuals))
outfile.close()

command = "%s/frappe1.1_linux64 frappeconfig.txt" % (frappe_dir)
os.system(command)
# args = shlex.split(command)
# p = os.popen(command,"r")
# while 1:
#     line = p.readline()
#     if not line: 
#     	break
# 	if flag:
# 		markers_to_remove.append()
# 	if line.startswith('The following problems have been identified'):
# 		flag=True

#     print line




# #open ped file and search for missing data
# for line in open(options.pedfile, 'r'):
# 	row = line.split(' ')
# 	count = 0
# 	count_missing = 0
# 	for item in row:

# 		if item == '0':
# 			print count
# 			count_missing +=1
# 			print 'achou '+str(count_missing)
			

# 		count +=1
		

# 	die()
	

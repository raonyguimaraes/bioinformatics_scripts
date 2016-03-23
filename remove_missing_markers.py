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
#open ped file and search for missing data
for line in open(options.pedfile, 'r'):
	row = line.split(' ')
	count = 0
	count_missing = 0
	for item in row:

		if item == '0':
			print count
			count_missing +=1
			print 'achou '+str(count_missing)
			

		count +=1
		

	die()
	

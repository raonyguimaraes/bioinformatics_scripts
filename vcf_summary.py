#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Description: VCF summary
from optparse import OptionParser
import os

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2011, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		

parser = OptionParser()

parser.add_option('-i', help='vcf', dest='vcffile') #, nargs=2
(options, args) = parser.parse_args()

vcffile=options.vcffile
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

readfile = open(vcffile, 'r')


### DATA ABOUT THE FILE
number_variant = 0
number_on_dbsnp = 0



for line in readfile:
    
    #header
    #if line.startswith('##'):
	#if line.startswith('##ApplyRecalibration'):
	    #print "Apply Recalibration Options\n\n\n"
	    
	    #recalibration_options = line.split('=', 1)[1].replace('"', '').split(' ')
	    #for option in recalibration_options:
		#print option
	#if line.startswith('##CombineVariants'):
	    
	    #combinevariants_options = line.split('=', 1)[1].replace('"', '').split(' ')
	    #for option in combinevariants_options:
		##print option
		#if option.startswith('variant'):
		    #print ''
		#elif option.startswith('rod_priority_list'):
		    #print ''
		#elif option.startswith('name'):
		    #print ''
		#elif option.startswith('source'):
		    #print ''
		#else:
		    #print option
	
	#if line.startswith('##FILTER'):
	    #print "\nFilter"
	    #filter_options = line.split('=', 1)[1][1:-2].split(',')
	    #for option in filter_options:
		#options = option.split('=', 1)
		#print "%s:%s" % (options[0], options[1])
	#if line.startswith('##FORMAT'):
	    #print "\n\nFormat"
	    #format_options = line.split('=', 1)[1][1:-2].split(',')
	    #for option in format_options:
		#options = option.split('=', 1)
		#print "%s:%s" % (options[0], options[1])
	#if line.startswith('##INFO'):
	    #print "\n\nINFO"
	    #format_options = line.split('=', 1)[1][1:-2].split(',')
	    #for option in format_options:
		#options = option.split('=', 1)
		#print "%s: " % (options[0]),
		#for item in options[1:]:
		    #print item,
		#print 
		
	    #die()
		
		
	    
	    #die()
	
    #only data
    if not line.startswith('#'):
	#count number of variants
	number_variant += 1
	variant = line.split('\t')
	#check if variant is at dbSNP
	if not variant[2].startswith('.'):
	    number_on_dbsnp += 1
	
	
	
	
	
	
percentage_in_dbsnp = round(float(number_on_dbsnp)/number_variant, 5)

print "Summary"
print "Number of variants: %s" % (number_variant)
print "Number of variants on dbSNP: %s %s %%" % (number_on_dbsnp, percentage_in_dbsnp)




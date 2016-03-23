#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Description: VCF summary
from optparse import OptionParser
import os


parser = OptionParser()

parser.add_option('-i', help='vcf', dest='vcffile') #, nargs=2

(options, args) = parser.parse_args()

vcffile=options.vcffile


filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

duplicates = 0
previous_position = 0
#For VCF One
vcf_one = open(vcffile, 'r')
for line in vcf_one:
    if not line.startswith('#'):
	line = line.split('\t')
	
	actual_position = line[1]
	if previous_position == actual_position:
	    print previous_line
	    print line
	    
	    duplicates += 1
	else:
	    previous_position = actual_position
	    previous_line = line
	    
print 'duplicates: %s' % (duplicates)
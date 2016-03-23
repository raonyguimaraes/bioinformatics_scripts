#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime



__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2011, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		

parser = OptionParser()

parser.add_option("-i", dest="vcf_file",
                  help="VCF File", metavar="VCF")
                  
(options, args) = parser.parse_args()

vcf_file = options.vcf_file
filename = vcf_file.split('/')[-1].split('.')[0]
print filename

chromossomes = []
for number in range(1,23):
    chromossomes.append(str(number))

chromossomes.append('X')
chromossomes.append('Y')

#open annovar file
vcf_file = open(vcf_file, 'r')

vcf_output =  open(filename+"_auto_sex.vcf","w")



for line in vcf_file:
    #print line,
    if line.startswith("#"):
	vcf_output.write(line)
    else:
	if line.split("\t")[0] in chromossomes:
	    vcf_output.write(line)


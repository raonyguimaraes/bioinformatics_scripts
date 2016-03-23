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
parser.add_option("-o", dest="out_vcf_file",
                  help="VCF File", metavar="VCF")                  
(options, args) = parser.parse_args()

vcf_file_path=options.vcf_file

#open annovar file
vcf_file = open(vcf_file_path, 'r')

vcf_output =  open(options.out_vcf_file,"w")


for line in vcf_file.readlines():
    #print line,
    if line.startswith("#"):
      #remove snpeff annotations
      if not line.startswith('##SnpEff'):
	vcf_output.write(line)
    else:
	vcf_output.write(line.replace('chr', ''))
	#if len(line.split("\t")[0]) <= 3:
	    #if line.startswith("MT") == False:
		
#command = "mv vcf_withoutchr.vcf %s" % (options.vcf_file)
#os.system(command)

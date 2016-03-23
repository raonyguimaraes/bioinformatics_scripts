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

vcf_file_path=options.vcf_file

#open annovar file
vcf_file = open(vcf_file_path, 'r')

filename = os.path.splitext(os.path.basename(str(options.vcf_file)))[0]

vcf_output =  open("%s.PASS.vcf" % (filename),"w")

for line in vcf_file:
	if line.startswith('#'):
		vcf_output.write(line)
	else:
		row = line.split('\t')
		if row[6] == ".":
			row[6] = "PASS"
		vcf_output.write("\t".join(row))
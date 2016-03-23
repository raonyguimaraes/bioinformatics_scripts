#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
import subprocess
import time


parser = OptionParser()

parser.add_option("-i", dest="vcffile",
                  help="VCF File", metavar="VCF")
parser.add_option("-o", dest="out_vcffile",
                  help="VCF File", metavar="VCF")
parser.add_option("-k", dest="k",
                  help="Number of SNPs", metavar="K")

(options, args) = parser.parse_args()
vcffile=options.vcffile
outfile = open(options.out_vcffile, 'w')

variant_list = []
for line in vcffile:
	if line.startswith('#'):
		outfile.writelines(line)
	else:
		variant_list.append(line)

while 
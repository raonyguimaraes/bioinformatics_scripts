#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
import csv

parser = OptionParser()
parser.add_option("-i", dest="vcf_file",
                  help="VCF File to Annotate", metavar="VCF")
                  
(options, args) = parser.parse_args()

output = open('sift_variants.txt', 'w')
for line in open(options.vcf_file):
    if not line.startswith('#'):
	line = line.split('\t')
	output.writelines('%s,%s,1,%s/%s\n' % (line[0], line[1], line[3], line[4]))
                  
#os.system("FILE1='%s'" % (''))
#os.system("FINALOUTDIR='%s'" % (os.getcwd()))
command = 'sh /lgc/programs/sift/SIFTwebservice/SIFTexome_nssnvs.sh /lgc/programs/sift/SIFTwebservice/SIFTexome_nssnvs.conf sift_variants.txt %s' % (os.getcwd())
os.system(command)
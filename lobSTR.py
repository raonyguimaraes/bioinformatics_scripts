#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime



__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, LobSTR"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Development"
		
#run example
#python lobSTR.py -i exome.sorted.bam

parser = OptionParser()
usage = "usage: %prog [options] -i exome.bam"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted in BAM format", metavar="BAM")
parser.add_option("-o", dest="output_file",
                  help="outputfile", metavar="out")
parser.add_option("-s", dest="sex",
                  help="M or F", metavar="SEX")
                  
                  
(options, args) = parser.parse_args()

input_file=options.input_file
output_file = options.output_file
sex = options.sex

lob_index = '/projects/rms_project/bin/lobstr_v1.0.6_2.linux64/bin'
index_lobSTR = '/lgc/datasets/lobSTR/index_trf_hg19/lobSTR_'

command = '%s/lobSTR -p 8 --bam -f %s --index-prefix %s -o %s' % (lob_index, input_file, index_lobSTR, output_file)
os.system(command)

#both
command = '%s/allelotype -v --command both --bam %s.aligned.bam --noise_model %s.noisemodel.txt --out %s --sex %s' % (lob_index, output_file, output_file, output_file, sex)
os.system(command)

#onlycall
command = '%s/allelotype -v --command simple --bam %s.aligned.bam --out %s --sex %s' % (lob_index, output_file, output_file, sex)
#os.system(command)



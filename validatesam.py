#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime



__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2013, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "2.3.4"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Develpment"
		
#run example
#python gatk.py -i alignment/exome.sorted.bam

parser = OptionParser()
usage = "usage: %prog [options] -i reads.bam"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted in BAM format", metavar="BAM")
# parser.add_option("-t", dest="target_array",
#                   help="Target Array", metavar="BEDFILE")
                  
(options, args) = parser.parse_args()

input_file=options.input_file
# target_array=options.target_array

pic_dir = '/lgc/programs/picard-tools-1.82'

command = 'java -jar %s/ValidateSamFile.jar I=%s  R=/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta VALIDATION_STRINGENCY=SILENT' % (pic_dir, input_file)
os.system(command)
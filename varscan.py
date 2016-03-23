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


parser = OptionParser()
usage = "usage: %prog [options] -i exome.bam"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted in BAM format", metavar="BAM")
parser.add_option("-o", dest="output_file",
                  help="outputfile", metavar="out")
                  
(options, args) = parser.parse_args()

input_file=options.input_file
output_file = options.output_file

reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.fasta"
varscan = '/lgc/programs/varscan/VarScan.v2.3.7.jar'
samtools_dir = '/lgc/programs/samtools-bcftools-htslib-1.0_x64-linux/bin'

command = '%s/samtools mpileup -f %s %s > %s.mpileup' % (samtools_dir, reference, input_file, input_file)
# os.system(command)

command = 'java -jar %s mpileup2indel %s.mpileup --output-vcf 1 > %s' % (varscan, input_file, output_file)
os.system(command)

# --min-coverage 30 \
# --min-reads2 10 \
# --min-avg-qual 30 \
# --p-value 0.99 \

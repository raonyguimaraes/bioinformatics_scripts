#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
		
#run example
#python bwa.py -f ../../input/exome/Ot729_index7_Exome_Sergio.Pena_RDNP_07262011-reads1-110716_I123_FCC043NABXX_L2_index7_1.fastq -q ../../input/exome/Ot729_index7_Exome_Sergio.Pena_RDNP_07262011-reads2-110716_I123_FCC043NABXX_L2_index7_2.fastq

parser = OptionParser()
usage = "usage: %prog [options] -f reads1.fastq -q reads2.fastq"
parser = OptionParser(usage=usage)

parser.add_option("-f", dest="reads1",
                  help="reads 1 in FASTQ format", metavar="FASTQ")
parser.add_option("-q", dest="reads2",
                  help="reads 2 in FASTQ format", metavar="FASTQ")
parser.add_option("-t", dest="target_array",
                  help="Target Array", metavar="BEDFILE")
                  
(options, args) = parser.parse_args()


reads1=options.reads1
reads2=options.reads2
target_array=options.target_array

print "Alignment"
command = 'python /lgc/scripts/bwa.py -f %s -q %s' % (reads1, reads2)
os.system(command)
print "GATK"
command = 'python /lgc/scripts/gatk.py -i exome.sorted.bam -t %s' % (target_array)
os.system(command)


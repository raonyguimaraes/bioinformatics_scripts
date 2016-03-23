#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime

#Example
#python /lgc/scripts/unifiergenotyper.py -i ../../data/Exome_11_AP.realigned-recalibrated.bam -t /lgc/datasets/exome_targets/Agilent_50MB/BED/029720_D_BED_20111101.bed


__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#run example
#python gatk.py -i alignment/exome.sorted.bam

parser = OptionParser()
usage = "usage: %prog [options] -f reads1.fastq -q reads2.fastq"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted in BAM format", metavar="BAM")
parser.add_option("-o", dest="output",
                  help="Output Filename", metavar="VCFFILE")
                  
(options, args) = parser.parse_args()

input_file=options.input_file

filename = ".".join(input_file.split("/")[-1].split(".")[:-1])


#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"
#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.test.fasta"
# reference="/lgc/datasets/gatk_data/2.5/human_g1k_v37.fasta"
#reference = '/lgc/datasets/gatk_data/b37/human_g1k_v37chr.fasta'
reference="/lgc/datasets/hg19/new/hg19/hg19.fa"
# reference = '/lgc/datasets/hg18/all/all.fasta'

pic_dir="/lgc/programs/picard-tools-1.128"

command = "java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/picard.jar ReorderSam I=%s O=karyotypic.%s REFERENCE=%s" % (pic_dir, input_file, input_file, reference)
# command =  """
# java -jar -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/AddOrReplaceReadGroups.jar I=%s O=rg%s LB=whatever PL=illumina PU=whatever SM=whatever
# """ % (pic_dir, input_file, input_file)      
os.system(command)    



      

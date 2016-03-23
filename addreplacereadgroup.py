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

reference="/lgc/datasets/hg19/new/hg19/hg19.fa"

dbsnp="/lgc/datasets/dbsnp/137/00-All.vcf"

pic_dir="/lgc/programs/picard-tools-1.109"
pic_dir="/storage2/programs/piccard_tools/picard-tools-1.119"

#GATK 1.5
# command = """
#       java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T UnifiedGenotyper \
#       -l INFO \
#       -I %s \
#       -R %s \
#       --dbsnp %s \
#       --genotype_likelihoods_model BOTH \
#       -stand_call_conf 50.0 \
#       -stand_emit_conf 0.0 \
#       -dcov 350 \
#       -A AlleleBalance \
#       -A DepthOfCoverage \
#       -A FisherStrand \
#       -o  %s.vcf\
#       -log %s-UnifiedGenotyper.log \
#       -nt 8
#       """ % (gatk_dir, input_file, reference, dbsnp, filename, filename)

command =  """
java -jar -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/AddOrReplaceReadGroups.jar I=%s O=rg%s LB=exome PL=illumina PU=exome SM=exome VALIDATION_STRINGENCY=SILENT
""" % (pic_dir, input_file, input_file)      
os.system(command)    



      

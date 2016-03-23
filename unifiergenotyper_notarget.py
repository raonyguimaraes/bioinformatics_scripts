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


reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"
reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.test.fasta"
# reference="/lgc/datasets/gatk_data/2.5/human_g1k_v37.fasta"
reference="/lgc/datasets/hg19/new/hg19/hg19.fa"


dbsnp="/projects/www/mendelmd_dev/annotator/data/dbsnp/00-All.vcf"

gatk_dir="/lgc/programs/GenomeAnalysisTK-1.5-21-g979a84a"
gatk_dir="/lgc/programs/GenomeAnalysisTK-3.1-1"

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

#GATK 2.6-4
command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T UnifiedGenotyper \
      -R %s \
      -I %s \
      -l INFO \
      -stand_call_conf 50.0 \
      -stand_emit_conf 10.0 \
      -dcov 200 \
      -A AlleleBalanceBySample \
      -o  %s \
      -log %s-UnifiedGenotyper.log \
      -nt 12
      """ % (gatk_dir, reference, input_file, options.output, options.output)
      #      --dbsnp %s \
# command = """cat %s.vcf | java -jar SnpSift.jar filter " ( QUAL >= 30 ) & (FILTER = 'PASS') & ( DP > 10 )" > filtered.vcf
# """ 
      # -A AlleleBalance \
      # -A DepthOfCoverage \
      # -A FisherStrand \
      
os.system(command)    



      

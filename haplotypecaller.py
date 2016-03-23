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
usage = "usage: %prog [options] -i exome1.bam exome2.bam"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted BAM format", metavar="BAM")
parser.add_option("-t", dest="target_array",
                  help="Target Array", metavar="BEDFILE")
parser.add_option("-o", dest="output",
                  help="Output Filename", metavar="VCFFILE")
                                    
                  
(options, args) = parser.parse_args()

input_file=options.input_file

filename = ".".join(input_file.split("/")[-1].split(".")[:-1])

target_array=options.target_array

#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"
#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.fasta"
#dbsnp="/lgc/datasets/dbsnp/dbsnp-135chr.vcf"

reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.test.fasta"
dbsnp="/lgc/datasets/dbsnp/141/All.vcf"

gatk_dir="/lgc/programs/gatk-protected/dist"
gatk_dir='/lgc/programs/GenomeAnalysisTK-3.3/'

#--dbsnp %s \

command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T HaplotypeCaller \
      -l INFO \
      -I %s \
      -R %s \
      --dbsnp %s \
      -stand_call_conf 50.0 \
      -stand_emit_conf 10.0 \
      -o  %s.vcf\
      -log %s-HaplotypeCaller.log \
      -L %s \
      """ % (gatk_dir, input_file, reference, dbsnp, options.output, filename, target_array)
#      -nct 4 \
      
os.system(command)    


#vcf_file = open("%s.vcf" % (filename), 'r')
#vcf_output = open(filename, 'w')
#for line in vcf_file :
    #if line.startswith("#"):
	#vcf_output.write(line)
    #else:
	#vcf_output.write(line.replace('chr', ''))
#vcf_output.close()
#vcf_file.close()

#command = 'mv %s %s.vcf' % (filename, filename)
#os.system(command)


      

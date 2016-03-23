#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime
import re
#Example
#python /projects/12exomes/data/unifiergenotyper.py -i /projects/12exomes/data/Exome_3_EDS.realigned-recalibrated.bam /projects/12exomes/data/Exome_4_ELS.realigned-recalibrated.bam /projects/12exomes/data/Exome_5_LS.realigned-recalibrated.bam /projects/12exomes/data/Exome_6_DC.realigned-recalibrated.bam -o quartet_only_SNPs_.vcf -t /lgc/datasets/exome_targets/SureSelect_All_Exon_V2_hg19.20110105.bed


__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#run example
#python unifiergenotyper.py -i /projects/12exomes/data/Exome_3_EDS.realigned-recalibrated.bam /projects/12exomes/data/Exome_4_ELS.realigned-recalibrated.bam /projects/12exomes/data/Exome_5_LS.realigned-recalibrated.bam /projects/12exomes/data/Exome_6_DC.realigned-recalibrated.bam

parser = OptionParser()
usage = "usage: %prog [options] -f reads1.fastq -q reads2.fastq"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="VCF File", metavar="VCF")
parser.add_option("-p", dest="ped_file",
                  help="Ped File", metavar="BEDFILE")
parser.add_option("-o", dest="output",
                  help="Output Filename", metavar="VCFFILE")
                  
(options, args) = parser.parse_args()

input_file=options.input_file



# filename = ".".join(input_file.split("/")[-1].split(".")[:-1])

#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"
#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.fasta"
#dbsnp="/lgc/datasets/dbsnp/dbsnp-135chr.vcf"

reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.fasta"
#dbsnp="/lgc/datasets/dbsnp/137/00-All.vcf"

gatk_dir="/lgc/programs/GenomeAnalysisTK-2.3-9-ge5ebf34"

#--dbsnp %s \

command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T PhaseByTransmission \
      -l INFO \
      -V %s \
      -R %s \
      -ped %s \
      -o  %s\
      -log %s-PhaseByTransmission.log \
      --pedigreeValidationType SILENT \
      """ % (gatk_dir, input_file, reference, options.ped_file, options.output, options.output)
#--output_mode EMIT_ALL_CONFIDENT_SITES \
print command
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


      
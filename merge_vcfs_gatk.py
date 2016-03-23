#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, Filter Analysis"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#run example
#python merge_vcfs_gatk.py -v vcf1 -q vcf2

parser = OptionParser()

parser.add_option("-v", dest="vcf_file_one",
                  help="VCF File", metavar="VCFa")
parser.add_option("-q", dest="vcf_file_two",
                  help="VCF File", metavar="VCFb")
parser.add_option("-o", dest="output_vcf",
                  help="VCF File", metavar="VCFc")

(options, args) = parser.parse_args()

gatk_dir='/lgc/programs/GenomeAnalysisTK-2.3-9-ge5ebf34'
reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"


#merge vcf files
command = '''java -Xmx40g -jar %s/GenomeAnalysisTK.jar \
   -R %s \
   -T CombineVariants \
   --variant:trio %s \
   --variant:1k %s \
   -o %s \
   --num_threads 8 \
   -genotypeMergeOptions UNIQUIFY''' % (gatk_dir, reference, options.vcf_file_one, options.vcf_file_two, options.output_vcf)

os.system(command)






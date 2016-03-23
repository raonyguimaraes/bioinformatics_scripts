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
#python gatk.py -i alignment/exome.sorted.bam

parser = OptionParser()

parser.add_option("-v", dest="vcf_file_one",
                  help="VCF File", metavar="VCFa")
parser.add_option("-q", dest="vcf_file_two",
                  help="VCF File", metavar="VCFb")
parser.add_option("-o", dest="output_vcf",
                  help="VCF File", metavar="VCFc")

(options, args) = parser.parse_args()

os.environ["PERL5LIB"] = "/lgc/programs/vcftools_0.1.10/lib/perl5/site_perl/"


#bgzip files
command = "bgzip %s" % (options.vcf_file_one)
os.system(command)

command = "bgzip %s" % (options.vcf_file_two)
os.system(command)

#tabix files
command = "tabix -p vcf %s.gz" % (options.vcf_file_one)
os.system(command)
command = "tabix -p vcf %s.gz" % (options.vcf_file_two)
os.system(command)

#merge-vcf
command = "/lgc/programs/vcftools_0.1.10/bin/vcf-merge %s.gz %s.gz > %s" % (options.vcf_file_one, options.vcf_file_two, options.output_vcf)
os.system(command)

#clean files
command = "rm -f %s.gz %s.gz.tbi" % (options.vcf_file_one, options.vcf_file_one)
os.system(command)
command = "rm -f %s.gz %s.gz.tbi" % (options.vcf_file_two, options.vcf_file_two)
os.system(command)





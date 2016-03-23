#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Description: VCF summary
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
		

parser = OptionParser()
#usage = python vcf_tools.py -i ~/projects/exome_analysis/input/boston/Boston_exome_Annotated.vcf
parser.add_option('-i', help='vcf', dest='vcffile') #, nargs=2
(options, args) = parser.parse_args()
  
vcffile=options.vcffile

filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

vcf_tools = '/lgc/programs/vcftools/bin'
command = "export PERL5LIB=/lgc/programs/vcftools/lib/perl5/site_perl/"
os.system(command)
#VCF-Validator
command = "%s/vcf-validator %s" % (vcf_tools, vcffile)
#os.system(command)

command = "%s/vcf-stats %s" % (vcf_tools, vcffile)
#os.system(command)

command = "%s/vcftools --vcf %s --freq --depth --het --hardy --TsTv-by-count --TsTv-by-qual --out %s " % (vcf_tools, vcffile, filename)
os.system(command)

#vcf-validator
#vcftools 
#vcf-stats
#!/usr/bin/env python
# -*- coding: utf-8 -*-
from optparse import OptionParser
import os
import csv

parser = OptionParser()
parser.add_option("-v", dest="vcf_file",
                  help="VCF File to Annotate", metavar="VCF")
#parser.add_option("-o", dest="vcf_output",
                  #help="VCF File to Annotate", metavar="VCF")

(options, args) = parser.parse_args()

vcffile = options.vcf_file

ind_file = open(vcffile)
for line in ind_file:
	if line.startswith('#CHROM'):
		row = line.strip().split('\t')
		individuals = row[9:]
		break
ind_file.close()

os.environ["PERL5LIB"] = "/lgc/programs/vcftools_0.1.10/lib/perl5/site_perl/"
for individual in individuals:
  print individual
  command = 'cat %s | /lgc/programs/vcftools_0.1.10/bin/vcf-subset -e -c %s | bgzip -c > %s.vcf.gz &' % (vcffile, individual, individual)
  os.system(command)
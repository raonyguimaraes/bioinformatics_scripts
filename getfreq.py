#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Description: VCF summary
from optparse import OptionParser
import os
import shlex, subprocess

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2011, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
frappe_dir = '/lgc/programs/frappe1.1/'

#1092exomes.sureselectv4.vcf.gz
import gzip


parser = OptionParser()

parser.add_option('-i', help='map', dest='inputfile') #, nargs=2

(options, args) = parser.parse_args()

inputfile=options.inputfile

filename = ".".join(inputfile.split("/")[-1].split(".")[:-1])

variant_list = []
for line in open(inputfile, 'r'):
	if line.startswith('#'):
		pass
	else:
		variant = line.split('\t')
	variant_list.append(variant[2])



print len(variant_list)
variants_dict = {}

#get freq for EUR, AMR,AFR
f = gzip.open('/projects/1000genomes/integrated_call_sets/1092exomes/1092exomes.EUR_AFR_AMR.sureselectv4.vcf.gz', 'rb')
for line in f:
	if line.startswith('#'):
		pass
	else:
		variant = line.split('\t')
		if variant[2] in variant_list:
			variant_id = variant[2]
			variants_dict[variant_id] = {}

			info = variant[7].split(';')
			for item in info:

				if item.startswith('AFR_AF'):
					variants_dict[variant_id]['AFR'] = item.split('=')[1]
				if item.startswith('AMR_AF'):
					variants_dict[variant_id]['AMR'] = item.split('=')[1]
				if item.startswith('EUR_AF'):
					variants_dict[variant_id]['EUR'] = item.split('=')[1]


#write file with freqs
output_file = open('allele_freq', 'w')
for variant in variants_dict:
	output_file.writelines('\t'.join(variants_dict[variant]['EUR'], variants_dict[variant]['AFR'], variants_dict[variant]['AMR'])+'\n')
	




file_content = f.read()
f.close()

#zcat aha!




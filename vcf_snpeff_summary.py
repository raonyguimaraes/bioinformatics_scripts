#!/usr/bin/env python
# -*- coding: utf-8 -*-
import scipy
import numpy as np
#Description: VCF summary
from optparse import OptionParser
import os
import rpy2.robjects as ro
from numpy import *
#import h5py
import csv

table = ro.r.table

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2011, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		

parser = OptionParser()

parser.add_option('-i', help='vcf', dest='vcffile') #, nargs=2
(options, args) = parser.parse_args()

vcffile=options.vcffile
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

readfile = open(vcffile, 'r')

#r_file = open('r_vcffile', 'w')

### DATA ABOUT THE FILE
number_variant = 0
number_on_dbsnp = 0
variants_list = []


def get_all_info_tags(vcffile):
    tempfile = open(vcffile, 'r')
    annotation_tags = set()
    for line in tempfile:
	if not line.startswith('#'):
	    variant = line.split('\t')
	    string = variant[7].split(';')
	    for element in string:
		element =  element.split('=')
		tag = element[0]	
		if tag not in annotation_tags:
		    annotation_tags.add(tag)
		    #annotation_tags.append(tag)
    tempfile.close()
    return annotation_tags

annotation_tags = get_all_info_tags(vcffile)

#r_file.writelines(", ".join(annotation_tags)+'\n')

print "Finished Looking for tags!!"


def parse_info_tag(string):
    string = string.split(';')
    information = {}
    for element in string:
	element =  element.split('=')
	tag = element[0]
	if len(element) > 1:
	    information[tag] = element[1]
	else:
	    information[tag] = 1
    #print information
    information_list = []
    for tag in annotation_tags:
	if tag in information:
	    information_list.append(str(information[tag]))
	else:
	    information_list.append('')
	    
    
    return information_list
    
for line in readfile:
    if not line.startswith('#'):
	#count number of variants
	number_variant += 1
	variant = line.strip().split('\t')
	#print len(variant[7].split(';'))
	information = parse_info_tag(variant[7])
	
	#variant = variant + information
	#print variant
	variant = variant[:7] + information + variant[8:]
	#print variant
	#die()
	
	
	#variant.append(information)
	#print len(variant)
	
	#r_file.writelines(", ".join(variant)+'\n')
	variants_list.append(variant)
	#np.append(variants_list, variant)
	#die()
	#check if variant is at dbSNP
	if not variant[2].startswith('.'):
	    number_on_dbsnp += 1
    if line.startswith('#CHROM'):
	vcf_header = line.strip().split('\t')

print "Finished Reading VCF"

variants_array = scipy.array(variants_list)
print "Array dimensions: rowsXcolumns"
print variants_array.shape
vcf_header = vcf_header[:7] + list(annotation_tags) + vcf_header[8:]

#print vcf_header
#write data to CSV File
c = csv.writer(open("variants_ramires.csv", "w"), quoting=csv.QUOTE_ALL)
c.writerow(vcf_header)
for row in variants_list:
  c.writerow(row)
c.close()



#vcf_header = vcf_header+list(annotation_tags)

#f = h5py.File('variants.h5', 'w')
#f['dset'] = variants_array
#f.close()



#Save parsed data for not having to process it anymore!!!!!




#index = vcf_header.index('SNPEFF_EFFECT')
#SNPEFF_EFFECT = variants_array[:, index]                       # get everything in the SNPEFF_EFFECT column

#array_test = SNPEFF_EFFECT
###########################start R analysis
#array_size = array_test.size
##convert row to column
#array_test = array_test.reshape(array_size,1)
#tlist = ro.StrVector(array_test)
#keyWordArgs = {'row.names':ro.StrVector(("seed"))}
#x = ro.r['as.data.frame'](table(tlist))
#ro.r['print'](x)

#DATA REFERENCE
#{'BaseQRankSum': '1.746', 'SNPEFF_FUNCTIONAL_CLASS': 'SILENT', 'SNPEFF_IMPACT': 'LOW', 'SNPEFF_TRANSCRIPT_ID': 'ENST00000377705', 'FS': '2.722', 'SNPEFF_GENE_NAME': 'NOL9', 'DP': '154', 'culprit': 'FS', 'SNPEFF_AMINO_ACID_CHANGE': 'P534', 'Dels': '0.00', 'HaplotypeScore': '4.7510', 'AC': '0', 'MQRankSum': '0.077', 'AF': '0.00', 'VQSLOD': '5.8077', 'AN': '2', 'MQ0': '0', 'SNPEFF_EXON_ID': 'exon_1_6592028_6592139', 'SNPEFF_GENE_BIOTYPE': 'protein_coding', 'ReadPosRankSum': '1.921', 'InbreedingCoeff': '-0.0085', 'HRun': '1', 'SNPEFF_EFFECT': 'SYNONYMOUS_CODING', 'MQ': '59.66', 'QD': '9.88', 'SNPEFF_CODON_CHANGE': 'ccG/ccA'}

#SNPEFF_FUNCTIONAL_CLASS


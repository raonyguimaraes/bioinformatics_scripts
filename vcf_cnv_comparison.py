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

parser.add_option('-v', help='vcf', dest='vcffile') #, nargs=2
parser.add_option('-b', help='vcf', dest='vcffile2') #, nargs=2
(options, args) = parser.parse_args()



(options, args) = parser.parse_args()

vcffile=options.vcffile
vcffile2=options.vcffile2

filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])
filename2 = ".".join(vcffile2.split("/")[-1].split(".")[:-1])

quality_threshold = 10.0


#For VCF One
vcf_one = open(vcffile, 'r')

vcf_one_snps = {}
### DATA ABOUT THE FILE
number_variants = 0
variants_dbsnp = 0
number_q30 = 0
print filename
for line in vcf_one:
    #only data
    if not line.startswith('#'):
		
	variant = line.strip().split('\t')
	
	snp = {}
	snp['raw'] = variant
	snp['id'] = '%s-%s' % (variant[0], variant[1])
	snp['chr'] = variant[0]
	snp['pos'] = variant[1]
	snp['snpid'] = variant[2]
	snp['qual'] = float(variant[5])
	snp['info'] = variant[-1]
	snp['cnv'] = 'chr%s:%s-%s' % (variant[0], variant[1], variant[-1].split(';')[-1].split('=')[-1].strip())
	
	if snp['qual'] > quality_threshold:
	    number_q30 += 1
	vcf_one_snps[snp['id']] = snp
	#count number of variants
	number_variants += 1
	
	



print "Summary"
print "Number of variants: %s" % (number_variants)
print 'Number of variants with QUAL>30: %s' % (number_q30)


print filename2
number_variants = 0
variants_dbsnp = 0
number_q30 = 0
vcf_two_snps = {}
vcf_two = open(vcffile2, 'r')
for line in vcf_two:
  
  if not line.startswith('#'):
		
	variant = line.strip().split('\t')
	
	
	snp = {}
	snp['raw'] = variant
	snp['id'] = '%s-%s' % (variant[0], variant[1])
	snp['chr'] = variant[0]
	snp['pos'] = variant[1]
	snp['snpid'] = variant[2]
	snp['qual'] = float(variant[5])
	snp['info'] = variant[-1]
	snp['cnv'] = 'chr%s:%s-%s' % (variant[0], variant[1], variant[-1].split(';')[-1].split('=')[-1].strip())
	#snp['genotype'] = variant[-1].split(':')[0]
	
	#count snps according to criteria
	#ignore 0/0 genotypes
	#if not snp['genotype'] == '0/0' or snp['genotype'] == '0|0':
	if snp['qual'] > quality_threshold:
	    number_q30 += 1
	vcf_two_snps[snp['id']] = snp
	#count number of variants
	number_variants += 1
	
	
	
	

print "Summary"
print "Number of variants: %s" % (number_variants)
print 'Number of variants with QUAL>30: %s' % (number_q30)


print 'VCF comparison'
#POSITIONS IN COMMON
genotypes_in_common = 0
genotypes_not_in_common = {}
for posid in vcf_one_snps:
    
    if posid in vcf_two_snps:
	
	if vcf_two_snps[posid]['info'] == vcf_one_snps[posid]['info']:
	    print 'boston, %s, %s' % (vcf_one_snps[posid]['cnv'], ", ".join(vcf_one_snps[posid]['raw']))
	    print 'ramires, %s, %s' % (vcf_two_snps[posid]['cnv'], ", ".join(vcf_two_snps[posid]['raw']))
	    
	    genotypes_in_common += 1
	else:
	    genotype = {}
	    genotype['vcfone'] = vcf_one_snps[posid]
	    genotype['vcftwo'] = vcf_two_snps[posid]
	    genotypes_not_in_common[posid] = genotype
	    
	
print 'Genotypes in common: %s' % (genotypes_in_common)

#for posid in genotypes_not_in_common:
    #print posid
    #print genotypes_not_in_common[posid]['vcfone']
    #print genotypes_not_in_common[posid]['vcftwo']

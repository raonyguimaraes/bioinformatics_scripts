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

quality_threshold = 0.0


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
		
	variant = line.split('\t')
	
	snp = {}
	snp['raw'] = variant
	snp['id'] = '%s-%s' % (variant[0].replace('chr', ''), variant[1])
	snp['chr'] = variant[0]
	snp['pos'] = variant[1]
	snp['snpid'] = variant[2]
	snp['ref'] = variant[3]
	snp['alt'] = variant[4]
	snp['qual'] = float(variant[5])
	snp['genotype'] = variant[-1].split(':')[0]
	#vcf_one_snps[snp['id']] = snp['genotype']
	#ignore 0/0 genotypes
	if not snp['genotype'] == '0/0' or snp['genotype'] == '0|0':
	    if snp['qual'] > quality_threshold:
		number_q30 += 1
		#vcf_one_snps[snp['id']] = snp['genotype']
		vcf_one_snps[snp['id']] = snp
		
	    #count number of variants
	    number_variants += 1
	    #check if variant is at dbSNP
	    if not variant[2].startswith('.'):
		variants_dbsnp += 1


percentage_in_dbsnp = round(float(variants_dbsnp)/number_variants, 5)
print "Summary"
print "Number of variants: %s" % (number_variants)
print "Number of variants on dbSNP: %s %s %%" % (variants_dbsnp, percentage_in_dbsnp)
print 'Number of variants with QUAL>30: %s' % (number_q30)


print filename2
number_variants = 0
variants_dbsnp = 0
number_q30 = 0
vcf_two_snps = {}
vcf_two = open(vcffile2, 'r')
for line in vcf_two:
  
  if not line.startswith('#'):
		
	variant = line.split('\t')
	
	
	snp = {}
	snp['raw'] = variant
	snp['id'] = '%s-%s' % (variant[0].replace('chr', ''), variant[1])
	snp['chr'] = variant[0]
	snp['pos'] = variant[1]
	snp['snpid'] = variant[2]
	snp['ref'] = variant[3]
	snp['alt'] = variant[4]
	snp['qual'] = float(variant[5])
	snp['genotype'] = variant[-1].split(':')[0]
	
	#count snps according to criteria
	#ignore 0/0 genotypes
	
	
	if not snp['genotype'] == '0/0' or snp['genotype'] == '0|0':
	    if snp['qual'] > quality_threshold:
		number_q30 += 1
		#vcf_two_snps[snp['id']] = snp['genotype']
		vcf_two_snps[snp['id']] = snp
	    #count number of variants
	    number_variants += 1
	    #check if variant is at dbSNP
	    if not variant[2].startswith('.'):
		variants_dbsnp += 1
	
	
	
percentage_in_dbsnp = round(float(variants_dbsnp)/number_variants, 5)
print "Summary"
print "Number of variants: %s" % (number_variants)
print "Number of variants on dbSNP: %s %s %%" % (variants_dbsnp, percentage_in_dbsnp)
print 'Number of variants with QUAL>30: %s' % (number_q30)

outputfile  = open('variants_in_common', 'w')

print 'VCF comparison'
genotypes_in_common = 0
genotypes_not_in_common = {}
for posid in vcf_two_snps:
    
    if posid in vcf_one_snps:
	#genotypes_in_common += 1
	if vcf_two_snps[posid]['genotype'] == vcf_one_snps[posid]['genotype']:
	    outputfile.write("\t".join([posid, vcf_one_snps[posid]['ref'], vcf_one_snps[posid]['alt'], vcf_one_snps[posid]['snpid'], vcf_two_snps[posid]['snpid'], vcf_one_snps[posid]['genotype'], vcf_two_snps[posid]['genotype']])+'\n')
	    
	    #outputfile.write("%s\t%s\t%s\n" % posid, vcf_one_snps[posid]['ref'], vcf_one_snps[posid]['alt'])
	    #print  vcf_two_snps[posid]
	    genotypes_in_common += 1
	    
	    
	    
	else:
	    genotype = {}
	    genotype['vcfone'] = vcf_one_snps[posid]
	    
	    genotype['vcftwo'] = vcf_two_snps[posid]
	    #print genotype
	    #die()
	    genotypes_not_in_common[posid] = genotype
	    
	
print 'Genotypes in common: %s' % (genotypes_in_common)

#for posid in genotypes_not_in_common:
    #print posid
    #print genotypes_not_in_common[posid]['vcfone']
    #print genotypes_not_in_common[posid]['vcftwo']

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
		
read_depth_treshold = 0
quality_treshold = 0.0
		
		
parser = OptionParser()

parser.add_option('-c', help='vcf', dest='child') #, nargs=2
parser.add_option('-p', help='vcf', dest='parents', nargs=2) #, nargs=2


(options, args) = parser.parse_args()

child=options.child

parents=options.parents

father_file = open(parents[0], 'r')
mother_file = open(parents[1], 'r')
child_file = open(child, 'r')
trio = {'father':{}, 'mother':{},'child':{}}

for line in father_file:
  if not line.startswith('#'):
    variant = line.split('\t')
    id = '%s-%s' % (variant[0].replace('chr', ''), variant[1])
    
    genotype_string = variant[-2].strip().split(':')
    if 'DP' in genotype_string:
      read_depth_index = genotype_string.index('DP')
      read_depth = int(variant[-1].strip().split(':')[read_depth_index])
      quality = float(variant[5])
      if read_depth >= read_depth_treshold:
	if quality >= quality_treshold:
	  genotype = variant[-1].strip().split(':')[0]
	  trio['father'][id] = genotype
print 'finished reading father'    
for line in mother_file:
  if not line.startswith('#'):
    variant = line.split('\t')
    id = '%s-%s' % (variant[0].replace('chr', ''), variant[1])
    genotype_string = variant[-2].split(':')
    if 'DP' in genotype_string:
      read_depth_index = genotype_string.index('DP')
      read_depth = int(variant[-1].strip().split(':')[read_depth_index])
      quality = float(variant[5])
      if read_depth >= read_depth_treshold:
	if quality >= quality_treshold:
	  genotype = variant[-1].strip().split(':')[0]
	  trio['mother'][id] = genotype
print 'finished reading mother'    
for line in child_file:
  if not line.startswith('#'):
    variant = line.split('\t')
    id = '%s-%s' % (variant[0].replace('chr', ''), variant[1])
    
    genotype_string = variant[-2].split(':')
    if 'DP' in genotype_string:
      read_depth_index = genotype_string.index('DP')
      read_depth = int(variant[-1].strip().split(':')[read_depth_index])
      quality = float(variant[5])
      if read_depth >= read_depth_treshold:
	if quality >= quality_treshold:
	  genotype = variant[-1].strip().split(':')[0]
	  trio['child'][id] = genotype
child_file.close()

print 'finished reading trio'
print len(trio['father'])
print len(trio['mother'])
print len(trio['child'])
    
child_denovo = {}    
#find denovo variants in child
for id in trio['child']:
  child_genotype = trio['child'][id]
  if id in trio['father']:
    if id in trio['mother']:
      #needs to be different from father
      if child_genotype != trio['father'][id] and child_genotype != trio['mother'][id]:
	#check for heterozygous
	if child_genotype == '0/0':
	  #check if both are not heterozygous
	  if trio['father'][id] != '0/1' and trio['mother'][id] != '0/1':
	    child_denovo[id] = trio['child'][id]	    
	elif child_genotype == '0/1':
	  #check if both are not homozygous
	  if trio['father'][id] != '1/1' and trio['mother'][id] != '1/1':
	    child_denovo[id] = trio['child'][id]
	  if trio['father'][id] != '0/0' and trio['mother'][id] != '0/0':
	    child_denovo[id] = trio['child'][id]
	else:
	  #1/1
	  #check if both are not heterozygous
	  if trio['father'][id] != '0/1' and trio['mother'][id] != '0/1':
	    child_denovo[id] = trio['child'][id]
	  
	#child_denovo[id] = trio['child'][id]

print 'denovo variants'	
print len(child_denovo)
#generate summary file
summary_file = open('denovo_variants_summary.txt', 'w')
header = 'variant_id\tchild\tfather\tmother\n'
summary_file.writelines(header)
for id in child_denovo:
  row = '%s\t%s\t%s\t%s\n' % (id, child_denovo[id], trio['father'][id], trio['mother'][id])
  summary_file.writelines(row)

denovo_vcf_file = open('denovo_variants.vcf', 'w')

child_file = open(child, 'r')
for line in child_file:
  if not line.startswith('#'):
    variant = line.split('\t')
    id = '%s-%s' % (variant[0].replace('chr', ''), variant[1])
    if id in child_denovo:
      denovo_vcf_file.writelines(line)
  else:
    denovo_vcf_file.writelines(line)
      
      
#generate vcf and summary file



  
  
  
  
  
  
    
    
  


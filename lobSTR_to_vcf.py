#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, Vaast Analysis"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#Options
parser = OptionParser()
parser.add_option('-i', help='tab', dest='tab')
(options, args) = parser.parse_args()

#1chromosome of matched STR
#2start position of matched STR
#3end position of matched STR
#4STR repeat
#5Period of STR repeat
#6Reference copy number (TRF notation)
#7Allelotype call (allele1, allele2) (given in # bp difference from reference) 
#8 Coverage
#9Number of reads supporting the reported genotype
#10Number of reads conflicting with the reported genotype
#11All alleles (given in # bp difference from reference:number of reads supporting that allele) that aligned to that locus separated by "/"
#12Score (ratio of likelihood of chosen allelotype to likelihood of all possibilities given the reads present). Note that the score is not always reliable and we currently do not recommend filtering on the score column. See "How should I filter my results to obtain only high quality genotypes" on the faq page.

print '''##fileformat=VCFv4.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta
##contig=<ID=1,length=249250621,assembly=b37>
##contig=<ID=10,length=135534747,assembly=b37>
##contig=<ID=11,length=135006516,assembly=b37>
##contig=<ID=12,length=133851895,assembly=b37>
##contig=<ID=13,length=115169878,assembly=b37>
##contig=<ID=14,length=107349540,assembly=b37>
##contig=<ID=15,length=102531392,assembly=b37>
##contig=<ID=16,length=90354753,assembly=b37>
##contig=<ID=17,length=81195210,assembly=b37>
##contig=<ID=18,length=78077248,assembly=b37>
##contig=<ID=19,length=59128983,assembly=b37>
##contig=<ID=2,length=243199373,assembly=b37>
##contig=<ID=20,length=63025520,assembly=b37>
##contig=<ID=21,length=48129895,assembly=b37>
##contig=<ID=22,length=51304566,assembly=b37>
##contig=<ID=3,length=198022430,assembly=b37>
##contig=<ID=4,length=191154276,assembly=b37>
##contig=<ID=5,length=180915260,assembly=b37>
##contig=<ID=6,length=171115067,assembly=b37>
##contig=<ID=7,length=159138663,assembly=b37>
##contig=<ID=8,length=146364022,assembly=b37>
##contig=<ID=9,length=141213431,assembly=b37>
##contig=<ID=X,length=155270560,assembly=b37>
##contig=<ID=Y,length=59373566,assembly=b37>
##INFO=<ID=PER,Number=1,Type=Integer,Description="Period of STR repeat">
##INFO=<ID=REF,Number=1,Type=Float,Description="Reference copy number (TRF notation)">
##INFO=<ID=ALL,Number=2,Type=String,Description="Allelotype call (allele1, allele2) (given in bp difference from reference)">
##INFO=<ID=COV,Number=1,Type=Integer,Description="Coverage">
##INFO=<ID=RSRG,Number=1,Type=Integer,Description="Number of reads supporting the reported genotype">
##INFO=<ID=RCRG,Number=1,Type=Integer,Description="Number of reads conflicting with the reported genotype">
##INFO=<ID=AA,Number=1,Type=String,Description="All alleles (given in # bp difference from reference:number of reads supporting that allele) that aligned to that locus separated by /">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample'''

filepath = open(options.tab, 'r')
for line in filepath:
  row = line.strip().split('\t')
  #print row
  info_list = []
  info_list.append('PER=%s' % row[4])
  info_list.append('REF=%s' % row[5])
  info_list.append('ALL=%s' % row[6])
  info_list.append('COV=%s' % row[7])
  info_list.append('RSRG=%s' % row[8])
  info_list.append('RCRG=%s' % row[9])
  info_list.append('AA=%s' % row[10])
  #genotype_info = '0/0:' % ()
  info_string = ";".join(info_list)
  print "\t".join([row[0].replace('chr', ''),row[1], '.', row[3], '.', row[-1], 'PASS', info_string, 'GT', '0/0'])
  
  #die()
  
  
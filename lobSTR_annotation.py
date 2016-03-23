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


filepath = open(options.tab, 'r')

#ccds = open('/lgc/datasets/CCDS.310112.txt', 'r')
#ccds_header = ccds.readline()
##print ccds_header
#gene_list = []
#for gene in ccds:
  #gene = gene.rstrip().split('\t')
  ##print gene
  #gene_obj = {}
  #gene_obj['name'] = gene[2]
  #gene_obj['chr'] = gene[0]
  #gene_obj['start'] = gene[7]
  #gene_obj['end'] = gene[8]
  #gene_obj['exons_list'] = gene[9][1:-1].split(',')
  ##print gene_obj['exons_list']
  #if gene_obj['exons_list'] != ['']:
    #gene_obj['exons'] = []
    #for exon in gene_obj['exons_list']:
      ##print exond
      #exon = exon.strip().split('-')
      #exon_obj = {}
      #exon_obj['start'] = exon[0]
      #exon_obj['end'] = exon[1]
      #gene_obj['exons'].append(exon_obj)
  
  #gene_list.append(gene_obj)
  ##gene_dict[gene[2]] = 
  
gene_list = []

#ref seq
#annotation = open('/lgc/datasets/ucsc-refseq.txt', 'r')
#knowngene
annotation = open('/lgc/datasets/ucsc-knowngene.tsv', 'r')

annotation_header = annotation.readline()

for gene in annotation:
  gene = gene.rstrip().split('\t')
  #print gene
  #die()
  #print gene
  gene_obj = {}
  gene_obj['name'] = gene[-4]
  gene_obj['chr'] = gene[2].replace('chr', '')
  
  gene_obj['start'] = gene[4]
  gene_obj['end'] = gene[5]
  
  gene_obj['exons_list_start'] = gene[9].split(',')
  gene_obj['exons_list_end'] = gene[10].split(',')
  
  #print gene_obj['exons_list']
  if gene_obj['exons_list_start'] != ['']:
    gene_obj['exons'] = []
    exon_count = 0
    for exon in gene_obj['exons_list_start']:
      
      #print exond
      exon = exon.strip().split('-')
      exon_obj = {}
      exon_obj['start'] = exon
      exon_obj['end'] = gene_obj['exons_list_start'][exon_count]
      
      gene_obj['exons'].append(exon_obj)
      
      exon_count += 1
  
  gene_list.append(gene_obj)
  #gene_dict[gene[2]] = 



#print len(gene_list)  
#die()
str_in_genes = 0
 
#for each STR
for line in filepath:
  row = line.strip().split('\t')
  
  chromossome = row[0].replace('chr', '')
  start = row[1]
  end = row[2]
  flag = False
  annotated = '\t'.join(row)
  #for each gene
  for gene in gene_list:
    #check chromossome
    if gene['chr'] == chromossome:
      #check if STR is in gene
      if start >= gene['start'] and end <= gene['end']:
	#for each exon
	
	flag = True
	annotated = annotated+'\t'+gene['name']
	str_in_genes += 1
	if 'exons' in gene:
	  count_exon = 0
	  for exon in gene['exons']:
	    count_exon += 1
	    if start >= exon['start'] and end <= exon['end']:
	      #print 'achou!'
	      #print gene, count_exon, row
	      
	      annotated = annotated+'\t'+str(count_exon)
	      
	      
  print annotated
  
          

print 'STRs in genes: %s' % (str_in_genes)

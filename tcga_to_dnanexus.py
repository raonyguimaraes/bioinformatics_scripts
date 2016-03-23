#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os


parser = OptionParser()
usage = "usage: %prog [options] -i file.tsv"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="TSV File ", metavar="TSV")
                  
(options, args) = parser.parse_args()

input_file=open(options.input_file, 'r')

filename = ".".join(options.input_file.split("/")[-1].split(".")[:-1])

outputfile = open("%s.dnanexus.csv" % filename, 'w')

header = input_file.readline().split('\t')
#print header


dnanexus_header = ['#var_index', 'chrom', 'left', 'right', 'ref_seq', 'var_type', 'zygosity', 'var_seq1', 'var_seq2', 'var_score', 'not_ref_score', 'coverage', 'read_count1', 'read_count2', 'conservation', 'gene_name', 'transcript_name', 'where_in_transcript', 'change_type1', 'ref_peptide1', 'var_peptide1', 'change_type2', 'ref_peptide2', 'var_peptide2', 'dbsnp', 'dbsnp_build']


tcag_header = ['chr', 'pos', 'ref', 'rsID', 'type', 'gt', 'typeOfMutation', '1000G allele frequency snp', '1000G ref|alt snp', '1000G allele frequency indel', '1000G ref|alt indel', 'complete genomics # of chromosomes called', 'complete genomics alt frequencies', 'complete genomics alt allele', 'TCAG exome # of chromosomes called', 'TCAG exome alt frequencies', 'TCAG exome alt allele', 'synonyms gene names', 'dbXrefs', 'gene description', 'RefSeq Effect impact', 'RefSeq Effect', 'RefSeq type of mutation', 'RefSeq codon change', 'RefSeq amino acid change', 'RefSeq gene name', 'RefSeq gene biotype', 'RefSeq coding', 'RefSeq transcript', 'RefSeq exon', '# of variants called in same gene', 'OMIM morbid map', 'OMIMgene map', 'NHLBI ESP EA MAF', 'NHLBI ESP AA MAF', 'NHLBI ESP ALL MAP', 'NHLBI ESP ALT Allele', 'A counts', 'C counts', 'G counts', 'T counts', 'N counts', 'DP', 'SB', 'QD', 'GATK filter', '# of potential het variants']


outputfile.write("\t".join(dnanexus_header)+'\n')

count = 0
for line in input_file:
    line = line.split('\t')
    output = []
    output.append(count)#var index
    count += 1
    
    output.append(line[0])#chrom
    output.append(line[1])#left
    output.append(int(line[1])+len(line[2])-1)#right
    output.append(line[2])#refseq
    #vartype
    ref = line[2]
    alt = line[5]
    vartype = line[4]
    mutation_type = line[6]
    if mutation_type == 'snp':
	output.append("SNP")#vartype
    else:#indel
	if len(alt) > len(ref):
	    output.append("INS")#vartype
	else:
	    output.append("DEL")#vartype
    
    output.append(line[4].title())#zygosity
    #var_seq
    if mutation_type == 'snp':
	output.append(line[5][0])#var_seq1
	output.append(line[5][1])#var_seq2
    else:
	if vartype == 'hom':
	    output.append(line[5])#var_seq1
	    output.append(line[5])#var_seq2
	else:
	    output.append(ref)#var_seq1
	    output.append(line[5])#var_seq2
    #varscore
    output.append(line[tcag_header.index('QD')])#varscore
    output.append(line[tcag_header.index('QD')])#notref_varscore
    
    output.append(line[tcag_header.index('DP')])#coverage
    try:
      output.append(line[tcag_header.index("%s counts" % line[5][0])])#read_count1
    except:
      output.append("")#read_count1
    try:
      output.append(line[tcag_header.index("%s counts" % line[5][1])])#read_count2
    except:
      output.append("")#read_count2
      
    output.append("")#conservation
    output.append(line[tcag_header.index("RefSeq gene name")])#gene_name
    output.append(line[tcag_header.index("RefSeq transcript")])#transcript_name
    
    output.append(line[tcag_header.index("RefSeq coding")].title())#where_in_transcript
    output.append(line[tcag_header.index("RefSeq Effect")].title())#change_type1
    

    output.append(line[tcag_header.index("RefSeq amino acid change")])#ref_peptide1
    output.append(line[tcag_header.index("RefSeq amino acid change")])#var_peptide1
    output.append(line[tcag_header.index("RefSeq Effect")].title())#change_type2
    output.append(line[tcag_header.index("RefSeq amino acid change")])#ref_peptide2
    output.append(line[tcag_header.index("RefSeq amino acid change")])#var_peptide2
    if line[3] == "0":
      output.append("")#dbsnp
    else:
      output.append(line[3])#dbsnp
    
    output.append("")#dbsnp_build
    new_output = []
    for item in output:
      new_output.append(str(item))
    output = new_output
      
    
    outputfile.write("\t".join(output)+'\n')

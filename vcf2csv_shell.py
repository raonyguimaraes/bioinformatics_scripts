#!/usr/bin/python
# -*- coding: utf-8 -*-

#script to convert vcf to csv
import os
import csv


from optparse import OptionParser

parser = OptionParser()
parser.add_option("-v", dest="vcf_file",
                  help="VCF File to Annotate", metavar="VCF")
(options, args) = parser.parse_args()

vcffile = options.vcf_file


#Get all annotation tags from a VCF File reading all file (lazy mode!)
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
    tempfile.close()
    annotation_tags = sorted(annotation_tags)
    
    return annotation_tags
    
    
def parse_info_tag(string, annotation_tags):
    string = string.split(';')
    information = {}
    for element in string:
	element =  element.split('=')
	tag = element[0]
	if len(element) > 1:
	    information[tag] = element[1]
	else:
	    information[tag] = 'Yes'
    information_list = []
    for tag in annotation_tags:
	if tag in information:
	    information_list.append(str(information[tag]))
	else:
	    information_list.append('')
    return information_list

def Get_vcfheader(filepath):
    vcffile = open(filepath, 'r')
    for line in vcffile:
	#print line
	if line.startswith("#CHROM"):
	    header_tags = line.strip().split('\t')
	if not line.startswith("#"):
	    break
    vcffile.close()
    return header_tags

vcf_header = Get_vcfheader(vcffile)

annotation_tags = get_all_info_tags(vcffile)

vcf_header = vcf_header[:7] +  annotation_tags + vcf_header[8:]

csvfilename = ".".join(os.path.basename(vcffile).split('.')[:-1])+'.csv'
csvfilepath = os.path.join(os.path.dirname(vcffile), csvfilename)

readfile = open(vcffile, 'r')
f = open(csvfilepath, "w")
csvfile = csv.writer(f, quoting=csv.QUOTE_ALL)
csvfile.writerow(vcf_header)
for line in readfile:
    if line.startswith('#CHROM'):
	vcf_header_original = line.strip().split('\t')
	vcf_header_original = vcf_header_original[:7] + list(annotation_tags) + vcf_header_original[8:]
    if not line.startswith('#'):
	variant = line.strip().split('\t')
	information = parse_info_tag(variant[7], annotation_tags)
	variant = variant[:7] + information + variant[8:]
	csvfile.writerow(variant)
	
	#new_variant = []
	##hack to old times
	
	#print vcf_header
	#print information
	#die()
	#for tag in vcf_header:
	    #tag_index = vcf_header_original.index(tag)
	    #new_variant.append(variant[tag_index])
	
f.close()

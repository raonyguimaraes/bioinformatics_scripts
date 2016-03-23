#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO

from optparse import OptionParser

parser = OptionParser()


parser.add_option("-i", dest="input_file",
                  help="Fastq in Illumina 1.3 Format", metavar="FASTQ")
                  
(options, args) = parser.parse_args()

input_file=options.input_file

#print input_file

in_handle = open(input_file)

#record_iterator = SeqIO.parse(fh, "fastq-illumina")
#out_handle = open("reads_python.fastq", "w")
#SeqIO.write(record_iterator, out_handle, "fastq")

out_handle = open("reads_python.fastq", "w")


count = SeqIO.convert(in_handle, "fastq-illumina", out_handle, "fastq")

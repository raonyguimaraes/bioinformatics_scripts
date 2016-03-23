#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime



__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2011, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		

parser = OptionParser()

parser.add_option("-i", dest="vcf_file",
                  help="VCF File", metavar="VCF")
parser.add_option("-o", dest="out_vcf_file",
                  help="VCF File", metavar="VCF")                  
(options, args) = parser.parse_args()

vcf_file=options.vcf_file

liftover_dir='/lgc/programs/gatk/public/perl/liftOverVCF.pl'
gatk_dir = ''

#sortbyref
#command = '/lgc/programs/gatk/public/perl/sortByRef.pl %s /lgc/datasets/hg18/all/all.fasta' % (vcf_file)
#sort vcf by hand



#Annotate with SNPEFF
filename = os.path.splitext(options.vcf_file)[0]

#Sort vcf file before annotation
# os.system("grep -E -v '^X|^Y|^MT|^#|^GL' %s.vcf | sort -n -k1 -k2 >> output.vcf" % (filename))
# os.system("grep '^#' %s.vcf > output.vcf" % (filename))
# os.system("grep -E '^X' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))
# os.system("grep -E '^Y' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))
# os.system("grep -E '^MT' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))

# chrM, chr1, chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13, chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, chrX, chrY,

os.system("grep '^#' %s.vcf > output.vcf" % (filename))


os.system("grep -E '^chrM' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))

for i in range(1,23):
  os.system("grep -E '^chr%s\t' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (i, filename))

os.system("grep -E '^chrX\t' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))
os.system("grep -E '^chrY' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))
# os.system("grep -E -v '^chrX|^chrY|^chrM|^#' %s.vcf | sort -n -k1 -k2 >> output.vcf" % (filename))

# os.system(command)

command = '%s -vcf output.vcf -gatk /lgc/programs/gatk/ -chain /lgc/datasets/hg18/hg18ToHg19.over.chain -newRef /lgc/datasets/gatk_data/hg19/ucsc.hg19 -oldRef /lgc/datasets/hg18/Homo_sapiens_assembly18 -out hg19%s.vcf' % (liftover_dir, filename)
print command

os.system(command)

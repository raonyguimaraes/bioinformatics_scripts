#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

parser.add_option('-i', help='vcf', dest='vcffile') #, nargs=2
(options, args) = parser.parse_args()

vcffile=options.vcffile
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

ann_dir="/lgc/programs/annovar"
#out_dir="output"

################################# Annovar ###############################
os.system('rm -rf annovar')
os.system('mkdir annovar')

command = "%s/convert2annovar.pl --format vcf4 --includeinfo %s > annovar/%s.annovar" % (ann_dir, vcffile, filename)
os.system(command)

#--ver1000g 1000g2011may
command = "%s/summarize_annovar.pl --ver1000g 1000g2012feb --buildver hg19 annovar/%s.annovar %s/humandb -outfile annovar/%s" % (ann_dir, filename, ann_dir, filename)
os.system(command)

#command = "%s/auto_annovar.pl --buildver hg19 --ver1000g 1000g2012feb -model dominant annovar/%s.annovar %s/humandb -outfile annovar/dominant.%s" % (ann_dir, filename , ann_dir, filename)
#os.system(command)


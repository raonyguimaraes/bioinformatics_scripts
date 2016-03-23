#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, Filter Analysis"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#run example
#python gatk.py -i alignment/exome.sorted.bam

parser = OptionParser()

parser.add_option("-p", dest="p1",
                  help="PED File", metavar="pedfile")
parser.add_option("-q", dest="p2",
                  help="PED File", metavar="pedfile")
parser.add_option("-o", dest="outfile",
                  help="PED File", metavar="pedfile")
(options, args) = parser.parse_args()

f1 = ".".join(options.p1.split("/")[-1].split(".")[:-1])
f1 = options.p1.replace('.ped','')

f2 = ".".join(options.p2.split("/")[-1].split(".")[:-1])
f2 = options.p2.replace('.ped','')
outfile = options.outfile

plink_dir = '/projects/relatedness/plink-1.07-x86_64'
#first identify the ones to remove
command = '%s/plink --file %s --merge %s.ped %s.map --recode --out %s --noweb --geno 0' % (plink_dir, f1, f2, f2, outfile)
os.system(command)

#commando remove snps
command = 'mv %s.missnp removesnps' % (outfile)
os.system(command)

print 'remove snps in file one'
command = '%s/plink --file %s --recode --out %s.snpsless --noweb --exclude removesnps' % (plink_dir, f1, f1)
os.system(command)

print 'remove snps in file two'
command = '%s/plink --file %s --recode --out %s.snpsless --noweb --exclude removesnps' % (plink_dir, f2, f2)
os.system(command)

print 'finally merge'
command = '%s/plink --file %s.snpsless --merge %s.snpsless.ped %s.snpsless.map --recode --out %s --noweb --geno 0' % (plink_dir, f1, f2, f2, options.outfile)
os.system(command)




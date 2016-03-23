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

parser.add_option('-i', help='vcf', dest='vcffile') #, nargs=2
(options, args) = parser.parse_args()

vcffile=options.vcffile
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

#vcf export
command = "export PERL5LIB=/lgc/programs/vcftools_0.1.10/lib/perl5/site_perl/"
os.system(command)

#calculate kinship
command = "/lgc/programs/vcftools_0.1.10/bin/vcftools --vcf %s --plink --out %s" % (vcffile, filename)
os.system(command)

command = "/projects/relatedness/plink-1.07-x86_64/plink --file %s --geno 0 --make-bed --noweb --out %s" % (filename, filename)
os.system(command)

command = "/projects/relatedness/plink/king -b %s.bed --kinship --prefix %s" % (filename,filename)
os.system(command)

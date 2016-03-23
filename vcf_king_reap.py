#!/usr/bin/env python
# -*- coding: utf-8 -*-

#Description: VCF summary
from optparse import OptionParser
import os
import shlex, subprocess

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
# control=options.control
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

vcf_dir = '/lgc/programs/vcftools_0.1.10/bin'

#remove chr
command = 'python /lgc/scripts/vcf_without_chr.py -i %s.vcf -o %s.nochr.vcf' % (filename, filename)
os.system(command)
#remove rs
command = 'python /lgc/scripts/remove_rs.py -i %s.nochr.vcf' % (filename)
os.system(command)
#remove y
command = 'python /lgc/scripts/remove_xy.py -i %s.nochr.nors.vcf' % (filename)
os.system(command)

#ped e map
command = "/lgc/programs/vcftools_0.1.10/bin/vcftools --vcf %s.nochr.nors.noxy.vcf --remove-indels --plink --out %s" % (filename, filename)
os.system(command)
die()

command = "/projects/relatedness/plink-1.07-x86_64/plink --file %s --geno 0 --make-bed --noweb --out %s" % (filename, filename)
os.system(command)
#king
# command = "/projects/relatedness/plink/king -b %s.bed --kinship --prefix %s" % (filename,filename)
command = 'python /lgc/scripts/remove_xy.py -i %s.nochr.nors.vcf' % (filename)
# os.system(command)

#reap
command = "python /lgc/scripts/reap_admixture_pedmap.py -p %s.ped -m %s.map" % (filename,filename)
os.system(command)


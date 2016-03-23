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

parser.add_option("-p", dest="ped",
                  help="PED File", metavar="pedfile")
parser.add_option("-m", dest="map",
                  help="MAP File", metavar="mapfile")

(options, args) = parser.parse_args()

pedfilename = ".".join(options.ped.split("/")[-1].split(".")[:-1])
mapfile = ".".join(options.map.split("/")[-1].split(".")[:-1])

plink_dir = '/projects/relatedness/plink-1.07-x86_64'

command = "/projects/relatedness/plink-1.07-x86_64/plink --file %s --maf 0.001 --geno 0 --output-missing-genotype 0 --make-bed --noweb --out %s" % (pedfilename, pedfilename)
os.system(command)

command = "/projects/relatedness/plink/king -b %s.bed --kinship --prefix %s" % (pedfilename,pedfilename)
os.system(command)

# command = "/projects/relatedness/plink/king -b %s.bed --kinship --ibs --related --prefix %s" % (pedfilename,pedfilename)
# os.system(command)
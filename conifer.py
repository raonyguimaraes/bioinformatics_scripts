#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime



__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, Conifer"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Development"
		
#run example
#python /lgc/scripts/conifer.py -i ../analysis270412/exome.sorted.dedup.real.fixed.recal.bam -t /lgc/programs/conifer_v0.2.1/probes.txt -o ramires


parser = OptionParser()
usage = "usage: %prog [options] -i exome.bam"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted in BAM format", metavar="BAM")
parser.add_option("-t", dest="target",
                  help="target", metavar="BED")
parser.add_option("-o", dest="output_file",
                  help="outputfile", metavar="out")
                  
(options, args) = parser.parse_args()

input_file=options.input_file
output_file = options.output_file
target = options.target

cnf_dir = '/lgc/programs/conifer_v0.2.2'

os.system('mkdir %s.RPKM' % output_file)
#Calculate RPKM from BAM
command = 'python %s/conifer.py rpkm --probes %s --input %s --output %s.RPKM/%s.rpkm.txt' % (cnf_dir, target, input_file, output_file, output_file)
os.system(command)

#cp rpkms from sample dir
os.system('cp %s/RPKM_data/* %s.RPKM' % (cnf_dir, output_file))

#Analyze
command = 'python %s/conifer.py analyze \
  	--probes %s \
  	--rpkm_dir %s.RPKM/ \
  	--output %s.analysis.hdf5 \
  	--svd 2 \
  	--write_svals %s.singular_values.txt \
  	--plot_scree screeplot.png \
  	--write_sd %s.sd_values.txt' % (cnf_dir, target, output_file, output_file, output_file, output_file)
os.system(command)

os.system('mkdir export_svdzrpkm')
#export all samples
command =  'python %s/conifer.py export \
  	--input %s.analysis.hdf5 \
  	--output ./export_svdzrpkm/' % (cnf_dir, output_file)
os.system(command)

#call
command = 'python %s/conifer.py call \
  	--input %s.analysis.hdf5 \
  	--output %s.calls.txt' % (cnf_dir, output_file,output_file)
os.system(command)

os.system('mkdir %s.call_imgs' % output_file)
#plot calls
command =  'python %s/conifer.py plotcalls \
  	--input %s.analysis.hdf5 \
  	--calls %s.calls.txt \
  	--outputdir ./%s.call_imgs/' % (cnf_dir, output_file, output_file, output_file)
os.system(command)
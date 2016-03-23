#!/usr/bin/env python
# -*- coding: utf-8 -*-

#run
#python /lgc/scripts/beagle.py -c ../analysis270412/exome.sorted.dedup.real.fixed.recal.raw.vcf -p ../../../12exomes/data/Exome_11_AP.var.annotated.vcf ../../../12exomes/data/Exome_12_IN.var.annotated.vcf

from optparse import OptionParser
import os

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, Beagle Analysis"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#Options
parser = OptionParser()
parser.add_option('-c', help='vcf', dest='child') #, nargs=2
parser.add_option('-p', help='vcf', dest='parents', nargs=2) #, nargs=2


(options, args) = parser.parse_args()

child=options.child
parents=options.parents

bg_dir = '/lgc/programs/beagle'
gatk_dir="/lgc/programs/GenomeAnalysisTK-1.5-21-g979a84a"
reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"
reference2='/lgc/datasets/hg19/all/all.fasta'

#without extension
child_filename = os.path.splitext(os.path.basename(child))[0]
parent0_filename = os.path.splitext(os.path.basename(parents[0]))[0]
parent1_filename = os.path.splitext(os.path.basename(parents[1]))[0]

#This walker is intended to be run after Beagle has successfully executed. The full calling sequence for using Beagle along with the GATK is:
command = 'java -Xmx2g -jar %s/GenomeAnalysisTK.jar -R %s -T ProduceBeagleInput -V %s -o %s' % (gatk_dir, reference, child, child_filename)
os.system(command)
command = 'java -Xmx2g -jar %s/GenomeAnalysisTK.jar -R %s -T ProduceBeagleInput -V %s -o %s' % (gatk_dir, reference2, parents[0], parent0_filename)
os.system(command)
command = 'java -Xmx2g -jar %s/GenomeAnalysisTK.jar -R %s -T ProduceBeagleInput -V %s -o %s' % (gatk_dir, reference2, parents[1], parent1_filename)
os.system(command)

#1. Run ProduceBeagleInputWalker.

#2. Run Beagle

#3. Uncompress output files

#4. Run BeagleOutputToVCFWalker.

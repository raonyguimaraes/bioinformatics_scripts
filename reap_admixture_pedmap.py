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


print 'merge with ped map from 1000genomes'
command = 'python /lgc/scripts/merge_pedmap.py -p %s -q /projects/1000genomes/integrated_call_sets/1092exomes/1092exomes_sureselect.ped -o %s_merged_ceu_yri_amr' % (pedfilename, pedfilename)
os.system(command)

print 'recode12'
filename2 = pedfilename+'_merged_ceu_yri_amr_12'
command = "/projects/relatedness/plink-1.07-x86_64/plink --file %s_merged_ceu_yri_amr --recode12 --output-missing-genotype 0 --missing --noweb --out %s" % (pedfilename, filename2)
os.system(command)

print 'tped'
command = '/projects/relatedness/plink-1.07-x86_64/plink --file %s --recode12 --output-missing-genotype 0 --transpose --out %s --noweb' % (filename2, filename2)
os.system(command)


# print command

print 'running admixture'

command = '/lgc/programs/admixture_linux-1.22/admixture %s.ped 3' % (filename2)
os.system(command)

#prep files for REAP
command = '''cut -d ' ' -f1,2 %s_merged_ceu_yri_amr.ped > myid.txt''' % (pedfilename)
os.system(command)

command = 'paste myid.txt %s.3.Q > admixturedata.proportions' % (filename2)
os.system(command)
command = '''sed 's/ /\t/g' admixturedata.proportions > adm.prop'''
os.system(command)

command = '/lgc/programs/reap/REAP/REAP -g %s.tped -p %s.tfam -a adm.prop -f %s.3.P -r 1 -k 3 -t -100' % (filename2, filename2,filename2)
os.system(command)


command = 'cat REAP_pairs_relatedness.txt | sort -r -n -k 9 > REAP_pairs_relatedness.ordered.txt'
os.system(command)


# command = "./REAP -g %s.tped -p %s.tfam -a trio1_bra_phasing_result.txt -f allelefreqle -k 2 -t 0.025 -r 2 -m" % (filename, filename)
# os.system(command)

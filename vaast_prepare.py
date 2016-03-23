#!/usr/bin/env python
# -*- coding: utf-8 -*-


#running program
#python /lgc/scripts/vaast_3targets_2parents.py -a vcf/ramires.snps.vcf vcf/boston.snps.vcf vcf/hungria.snps.vcf -p vcf/Exome_11.snps.vcf vcf/Exome_12.snps.vcf
#rm -rf cdr gvf simulations vat


from optparse import OptionParser
import os

__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, Vaast Analysis"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#Options
parser = OptionParser()
parser.add_option('-i', help='vcf', dest='individual') #, nargs=2

(options, args) = parser.parse_args()

individual=options.individual


#VAAST ANALYSIS 1.0.1
#vaast_dir="/projects/trio_analysis/bin/VAAST_Code_1.0.1/bin"
vaast_dir="/lgc/programs/VAAST_Code_1.0.3/bin"

out_dir="output"

genes="/lgc/datasets/vaast_data2/hg19/Features/refGene_hg19.gff3"
reference="/lgc/datasets/vaast_data2/hg19/Fasta/chrAll.fasta"
background="/lgc/datasets/vaast_data2/hg19/Background_CDR/1KG_179-refGene-hg19-liftover.cdr"
#background="/lgc/datasets/vaast_data2/hg19/Background_CDR/189individuals.cdr"


os.system("mkdir gvf")
os.system("mkdir simulations")
os.system("mkdir vat")
os.system("mkdir cdr")
#Convert vcf files to GVF Files
#command = "perl %s/vaast_tools/vaast_converter -f VCF --path gvf/ -b hg19 %s %s" % (vaast_dir, " ".join(affecteds), " ".join(parents))

command = "perl %s/vaast_tools/vaast_converter -f VCF --path gvf/ -b hg19 %s" % (vaast_dir, individual) #
os.system(command)
individuals = []
path = os.getcwd()+"/gvf"
for fname in os.listdir(path):
    if fname.endswith('.gvf'):
      individual = {}
      individual['name'] = fname.split('.')[0]
      individual['path'] = path+'/'+fname
      individuals.append(individual)


##VAT
print "Start VAT"
command = "%s/VAT --gender male -f %s -a %s -b hg19 %s > vat/%s.vat.gvf" % (vaast_dir, genes, reference, individuals[0]['path'], individuals[0]['name'])
os.system(command)
        



print "Start VST"

#VST Ramires
command = "%s/VST -o 'U(0)' -b hg19 vat/%s.vat.gvf > cdr/individual.cdr" % (vaast_dir, individuals[0]['name'])
os.system(command)
    
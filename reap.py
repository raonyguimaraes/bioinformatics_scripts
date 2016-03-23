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
		
frappe_dir = '/lgc/programs/frappe1.1/'

#1092exomes.sureselectv4.vcf.gz


parser = OptionParser()

parser.add_option('-i', help='vcf', dest='vcffile') #, nargs=2
parser.add_option('-d', help='gzvcf', dest='control')

(options, args) = parser.parse_args()

vcffile=options.vcffile
control=options.control
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

vcf_dir = '/lgc/programs/vcftools_0.1.10/bin'

#vcf export
command = "export PERL5LIB=/lgc/programs/vcftools_0.1.10/lib/perl5/site_perl"
os.system(command)

#compress
command = '/lgc/programs/tabix-0.2.6/bgzip %s' % (vcffile)
# os.system(command)
#index
command = '/lgc/programs/tabix-0.2.6/tabix -p vcf %s.gz' % (vcffile)
# os.system(command)


#merge vcf from input with 1000genomes
command = '%s/vcf-merge %s.vcf.gz %s | bgzip -c > %s.ceu_afr_amr.vcf.gz' % (vcf_dir, filename, control, filename)
# os.system(command)
#index
command = 'tabix -p vcf %s.ceu_afr_amr.vcf.gz' % (filename)
# os.system(command)

# die()
#remove missing genotypes
command = '''/lgc/programs/vcftools_0.1.10/bin/vcftools 
		 --remove-indels \
		 --gzvcf %s.ceu_afr_amr.vcf.gz \
         --max-alleles 2 \
         --min-alleles 2 \
         --geno 1 \
         --IMPUTE \
         --phased \
         --remove-indels \
         --recode --out  %s.ceu_afr_amr.no_missing.vcf''' % (filename, filename)
# print command
# os.system(command)

# die()
#generate tped
command = "/lgc/programs/vcftools_0.1.10/bin/vcftools --vcf %s.ceu_afr_amr.no_missing.vcf.recode.vcf --plink-tped --out %s.eur_afr_amr" % (filename, filename)
# os.system(command)

#generate ped
command = "/lgc/programs/vcftools_0.1.10/bin/vcftools --vcf %s.ceu_afr_amr.no_missing.vcf.recode.vcf --plink --out %s.eur_afr_amr" % (filename, filename)
# os.system(command)

#recode
filename2 = filename+'_12'
command = "/projects/relatedness/plink-1.07-x86_64/plink --file %s.eur_afr_amr --recode12 --missing --noweb --out %s" % (filename, filename2)
# os.system(command)
# print command

# command = "/projects/relatedness/plink-1.07-x86_64/plink --file %s --recode --missing --geno 0.1 --exclude remove_in_chry.txt --noweb --out %s" % (filename2, filename2)
# os.system(command)
output = subprocess.check_output("wc -l %s.ped" % (filename2), shell=True)
n_individuals = int(output.split()[0])
output = subprocess.check_output("wc -l %s.map" % (filename2), shell=True)
n_markers = int(output.split()[0])

#frappe
outfile = open("frappeconfig.txt", 'w')
outfile.writelines("""GenotypeFile  ="%s.ped" ## Mandatory genotype file name
MaxIter=  10000    ## maximum number of EM iterations
K      =  3  ## Number of ancestral individuals
M      = %s  ## Number of markers
I      =  %s  ## Number of individuals;
step   =  200 ### Define either step or Nout
Nout   =  0
printP  = 1     ## optional
threshold= 10000  ## optional convergence threshold

""" % (filename2,n_markers,n_individuals))
outfile.close()

command = "%s/frappe1.1_linux64 frappeconfig.txt" % (frappe_dir)
os.system(command)

# args = shlex.split(command)
# p = os.popen(command,"r")
# while 1:
#     line = p.readline()
#     if not line: 
#     	break
# 	if flag:
# 		markers_to_remove.append()
# 	if line.startswith('The following problems have been identified'):
# 		flag=True

#     print line


# p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# p = subprocess.call(args, stdout=subprocess.STDOUT,shell=True)
# output, errors = p.communicate()
# print output, errors
# print p
# die()

# for line in output:
# 	print 'hello'
# 	print line




# genomes1k_file = open()
#create allele frequency file
# all_freq_file = open('allelefreqle.txt', 'w')



# command = "/projects/relatedness/plink/king -b %s.bed --kinship --prefix %s" % (filename,filename)
# os.system(command)



command = "./REAP -g %s.tped -p %s.tfam -a trio1_bra_phasing_result.txt -f allelefreqle -k 2 -t 0.025 -r 2 -m" % (filename, filename)
os.system(command)

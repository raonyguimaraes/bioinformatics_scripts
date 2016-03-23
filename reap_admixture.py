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
# parser.add_option('-d', help='gzvcf', dest='control')


(options, args) = parser.parse_args()

vcffile=options.vcffile
# control=options.control
control = ' /projects/1000genomes/integrated_call_sets/1092exomes/1092exomes.EUR_AFR_AMR.sureselectv4.vcf.gz'
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

vcf_dir = '/lgc/programs/vcftools_0.1.10/bin'

#vcf export
command = "export PERL5LIB=/lgc/programs/vcftools_0.1.10/lib/perl5/site_perl"
os.system(command)

print 'filter by QUAL, PASS, READDepth'
command = '%s/vcftools --vcf %s --out %s.snps.q30.rd30 --recode --remove-indels' % (vcf_dir, vcffile, filename) #--minQ 30 --min-meanDP 30.0 1000exomes
os.system(command)

infile = open('%s.snps.q30.rd30.recode.vcf' % (filename), 'r')
outfile = open('%s.snps.q30.rd30.pass.vcf' % (filename), 'w')
#filter by pass
for line in infile:
   if line.startswith('#'):
      outfile.writelines(line)
   else:
      variant = line.split('\t')
      if variant[0] in ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']:
         if variant[6] in ['.', 'PASS']:
            outfile.writelines(line)


infile.close()
outfile.close()

filtered_file = '%s.snps.q30.rd30.pass.vcf' % (filename)

print 'merging vcfs'
#merge with gatk
command = '''python /lgc/scripts/merge_vcfs_gatk.py \
            -v %s.snps.q30.rd30.pass.vcf \
            -q /projects/1000genomes/integrated_call_sets/1092exomes/1092exomes.EUR_AFR_AMR.sureselectv4.vcf \
            -o %s.merged.with_ceu_afr_amr.vcf''' % (filename, filename)

os.system(command)

#compress
command = '/lgc/programs/tabix-0.2.6/bgzip %s.merged.with_ceu_afr_amr.vcf' % (filename)
os.system(command)

#index
command = 'tabix -p vcf %s.merged.with_ceu_afr_amr.vcf.gz' % (filename)
os.system(command)

print 'remove missing genotypes'
command = '''/lgc/programs/vcftools_0.1.10/bin/vcftools \
		 --gzvcf %s.merged.with_ceu_afr_amr.vcf.gz \
         --max-alleles 2 \
         --min-alleles 2 \
         --geno 1 \
         --phased \
         --remove-indels \
         --recode --out  %s.ceu_afr_amr.no_missing''' % (filename, filename)
# print command
os.system(command)

#generate ped
command = "/lgc/programs/vcftools_0.1.10/bin/vcftools --vcf %s.ceu_afr_amr.no_missing.recode.vcf --plink --out %s.eur_afr_amr" % (filename, filename)
os.system(command)

#recode
filename2 = filename+'_12'
command = "/projects/relatedness/plink-1.07-x86_64/plink --file %s.eur_afr_amr --recode12 --output-missing-genotype 0 --missing --noweb --out %s" % (filename, filename2)
os.system(command)

#tped
command = '/projects/relatedness/plink-1.07-x86_64/plink --file %s --recode12 --output-missing-genotype 0 --transpose --out %s --noweb' % (filename2, filename2)
os.system(command)


# print command

print 'running admixture'

command = '/lgc/programs/admixture_linux-1.22/admixture %s.ped 3' % (filename2)
os.system(command)

#prep files for REAP
command = 'cut -f1,2 %s.eur_afr_amr.ped > myid.txt' % (filename)
os.system(command)

command = 'paste myid.txt %s.3.Q > admixturedata.proportions' % (filename2)
os.system(command)
command = '''sed 's/ /\t/g' admixturedata.proportions > adm.prop'''
os.system(command)

command = '/lgc/programs/reap/REAP/REAP -g %s.tped -p %s.tfam -a adm.prop -f %s.3.P -r 1 -k 3' % (filename2, filename2,filename2)
os.system(command)


command = 'cat REAP_pairs_relatedness.txt | sort -r -n -k 9 > REAP_pairs_relatedness.ordered.txt'
os.system(command)


# command = "./REAP -g %s.tped -p %s.tfam -a trio1_bra_phasing_result.txt -f allelefreqle -k 2 -t 0.025 -r 2 -m" % (filename, filename)
# os.system(command)

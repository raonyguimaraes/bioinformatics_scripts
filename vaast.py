#!/usr/bin/env python
# -*- coding: utf-8 -*-


#running program
#python vaast.py -d ../yoruban_analysis/NA19240.YRI.hg19.with_sickle_cell_disease.vcf -f ../yoruban_analysis/NA19239.YRI.hg19.with_sickle_cell_disease.vcf -m ../yoruban_analysis/NA19238.YRI.hg19.with_sickle_cell_disease.vcf -a ../yoruban_analysis/NA18486.with_sickle_cell_disease.YRI.hg19.vcf

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
		
#Options
parser = OptionParser()
parser.add_option('-i', help='vcf', dest='affected') #, nargs=2

(options, args) = parser.parse_args()

affected=options.affected

#start filenames always with "NA76298" tag

affected_filename = ".".join(affected.split("/")[-1].split(".")[:-1])
#print affected_filename
#die()

out_dir="output"

#VAAST ANALYSIS 1.0.1
vaast_dir="/projects/trio_analysis/bin/VAAST_Code_1.0.1/bin"
out_dir="output"

genes="/projects/trio_analysis/input/vaast_data/hg19/Features/refGene_hg19.gff3"
reference="/projects/trio_analysis/input/vaast_data/hg19/Fasta/chrAll.fasta"
background="/projects/trio_analysis/input/vaast_data/hg19/Background_CDR/1KG_179-refGene-hg19-liftover.cdr"

#GENERATE A Background with 189 Genomes and with 1029 genomes


#select only variants from Chr1-22,X,Y
#command "~/lgc/scripts/vcf_withchr.py -i ../gatk/exome.sorted.dedup.real.fixed.recal.raw.vcf"
#os.system(command)

#Convert vcf files to GVF Files
command = "perl %s/vaast_tools/vaast_converter -b hg19 %s" % (vaast_dir, affected)
#os.system(command)

##VAT
print "Start VAT"
#affected
command = "%s/VAT -f %s -a %s -g male -b hg19 Exome.gvf > %s.vat.gvf" % (vaast_dir, genes, reference, affected_filename)
#os.system(command)


print "Start VST"

##VST Daugther
command = "%s/VST -b hg19 -o 'U(0)' %s.vat.gvf > %s.cdr" % (vaast_dir, affected_filename, affected_filename)
#os.system(command)

############## Finally Running VAAST !!!


print "Start VAAST"

#command = "time %s/VAAST -iht r  -lh n  -fast_gp  -d 1e6 -o  %s -m lrt  -k %s %s %s.cdr" % (vaast_dir, affected_filename, genes, background, affected_filename)
command = "time %s/VAAST -mp1 5 -mp2 5  -iht r  -lh y -pnt c  -fast_gp  -m lrt -d 1e6 -k --outfile %s %s %s Exome.cdr" % (vaast_dir, affected_filename, genes, background)
#print command
os.system(command)

command = "~/projects/trio_analysis/bin/VAAST_Code_1.0.1/bin/vaast_tools/simple_output.pl vcf_withchr.vaast > vaast.genelist"
os.system(command)

#INFO : processing_finished : Output can be found in vcf_withchr.vaast
#106844.33user 3080.28system 5:22:43elapsed 567%CPU (0avgtext+0avgdata 131395648maxresident)k
#0inputs+9273720outputs (0major+13023585minor)pagefaults 0swaps




#command = "%s/VAAST -mp1 8 --mode lrt --outfile %s %s %s %s.cdr" % (vaast_dir, affected_filename, genes, background, affected_filename)
#os.system(command)

##affected
#command = "time %s/VAAST -mp1 8 -mp2 8 --mode lrt --codon_bias --gp 10000 --outfile %s/vaast_results/affected_analysis %s %s %s/affecteds.cdr 2>run2_log" % (vaast_dir, out_dir, genes, background, out_dir)
##os.system(command)



##Parametrization
#command = "%s/VAAST -t %s/parents_union.cdr -pnt c -iht r -lh n -fast_gp -d 1e10 --significance 2.4e-6 -mp1 4 -mp2 4 -o trio_nalaysis_v1 -r 0.00035 -m lrt -k %s %s %s/daughter.cdr" % (vaast_dir, out_dir, genes, background, out_dir)
##os.system(command)
##generate list
#command = "%s/vaast_tools/simple_output.pl %s/trio_nalaysis_v1 > trio_nalaysis_v1.genelist" % (vaast_dir, out_dir)
##os.system(command)

#command = "time %s/VAAST -t %s/parents_union.cdr -pnt i -iht r -lh n -fast_gp -d 1e10 --significance 2.4e-6 -mp1 4 -o trio_nalaysis_v1_pnt_i -r 0.00035 -m lrt -k %s %s %s/daughter.cdr" % (vaast_dir, out_dir, genes, background, out_dir)
##os.system(command)
##generate list
#command = "%s/vaast_tools/simple_output.pl %s/trio_nalaysis_v1_pnt_i > trio_nalaysis_v1_pnt_i.genelist" % (vaast_dir, out_dir)
##os.system(command)



#command = "time %s/VAAST --mode lrt --codon_bias --inheritance r --penetrance i --gp 10000 --trio %s/parents_union.cdr --outfile %s/trio_analysis_with_179bg_withpenetrance %s %s %s/daughter.cdr" % (vaast_dir, out_dir, out_dir, genes, background, out_dir)
##os.system(command)
##generate list
#command = "%s/vaast_tools/simple_output.pl %s/trio_analysis_with_179bg_withpenetrance > trio_analysis_with_179bg_withpenetrance.genelist" % (vaast_dir, out_dir)
##os.system(command)









#%s/daughter.cdr







#VAAST
#time $VAAST_DIR/VAAST -p 8 -q 8 --mode lrt --rate 0.001 --codon_bias --gp 10000 --outfile $OUT_DIR/normal $GENES $BACKGROUND $OUT_DIR/affecteds.cdr

#Trio Analysis

#Generate parents.cdr using intersection or union ? (Intersection)


#time $VAAST_DIR/VAAST --mode lrt --codon_bias --inheritance r --penetrance i --gp 10000 --trio $OUT_DIR/parents.cdr --outfile $OUT_DIR/vaast_output/trio_analysis_with_179bg_withpenetrance $GENES $BACKGROUND $OUT_DIR/daughter.cdr

#time $VAAST_DIR/VAAST --mode lrt --codon_bias --gp 10000 --trio $OUT_DIR/parents.cdr --outfile $OUT_DIR/trio_analysis $GENES $BACKGROUND $OUT_DIR/daughter.cdr

#VAAST -m lrt –iht r –lh n –pnt c -r 0.00035 -p 4 -w 100000 -d 100000 -k -g 4 -x 35bp_se.wig -o <output_file> <feature_file> <189genome.cdr> <target> 

#AAS with OMIM command-line: 

#VAAST -q <degree of parallelization> -f -c <chromosome_number> -g 4 -k -d 100000000 -j 2.38E-6 -m lrt -o <output_file> <feature_file > <189genome.cdr> <target> 

#AAS with blosum62 command-line: 

#VAAST -q <degree of parallelization> -f -c <chromosome_number> -g 4 -k -b blosum62.matrix -d 100000000 -j 2.38E-6 -m lrt -o <output_file> <feature_file> <189genome.cdr> <target> 

#No AAS command-line: 

#VAAST -q <degree of parallelization> -f -c <chromosome_number> -g 4 -d 100000000 -j 2.38E-6 -m lrt -o <output_file> <feature_file> <189genome.cdr> <target> 

#WSS command-line: 

#VAAST -q <degree of parallelization> -c <chromosome_number> -d 100000000 -j 2.38E-6 -m wss -o <output_file> <feature_file> <189genome.cdr> <target>

#$VAAST_DIR/VAAST --mode lrt --codon_bias --gp 10000 --outfile $OUT_DIR/2affected $GENES $BACKGROUND $OUT_DIR/affecteds.cdr




#try this
#VAAST -m lrt –iht r –lh n –pnt c -r 0.00035 -p 4 -w 100000 -d 100000 -k -g 4 -x 35bp_se.wig -o <output_file> <feature_file> <189genome.cdr> <target> 

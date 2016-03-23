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
parser.add_option('-a', help='vcf', dest='affecteds', nargs=3) #, nargs=2
parser.add_option('-p', help='vcf', dest='parents', nargs=2) #, nargs=2


(options, args) = parser.parse_args()

affecteds=options.affecteds
parents=options.parents

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

for individual in affecteds:
  command = "perl %s/vaast_tools/vaast_converter -f VCF --path gvf/ -b hg19 %s" % (vaast_dir, individual) #
  #os.system(command)

for individual in parents:
  command = "perl %s/vaast_tools/vaast_converter -f VCF --path gvf/ -b hg19 %s" % (vaast_dir, individual)#
  #os.system(command)

#die()

#get files in GVF Folter
individuals = []
path = os.getcwd()+"/gvf"
for fname in os.listdir(path):
    if fname.endswith('.gvf'):
      individual = {}
      individual['name'] = fname.split('.')[0]
      individual['path'] = path+'/'+fname
      individuals.append(individual)
    
#print individuals
#die()
##VAT
print "Start VAT"
#affecteds
for individual in individuals:
    if individual['name'] == '0277-P1_1' or individual['name'] == 'Exome_12_IN':
      command = "%s/VAT --gender female -f %s -a %s -b hg19 %s > vat/%s.vat.gvf" % (vaast_dir, genes, reference, individual['path'], individual['name'])
    else:
      command = "%s/VAT --gender male -f %s -a %s -b hg19 %s > vat/%s.vat.gvf" % (vaast_dir, genes, reference, individual['path'], individual['name'])
    #print command
    #os.system(command)
    



print "Start VST"

#VST Ramires
command = "%s/VST -o 'U(0)' -b hg19 vat/%s.vat.gvf > cdr/ramires.cdr" % (vaast_dir, individuals[0]['name'])
#os.system(command)

##VST Affecteds
command = "%s/VST -o 'U(0..2)' -b hg19 vat/%s.vat.gvf vat/%s.vat.gvf vat/%s.vat.gvf > cdr/affecteds.cdr" % (vaast_dir, individuals[0]['name'], individuals[2]['name'], individuals[4]['name'])
#os.system(command)
##VST Parents
command = "%s/VST -o 'U(0..1)' -b hg19 vat/%s.vat.gvf vat/%s.vat.gvf > cdr/parents.cdr" % (vaast_dir, individuals[1]['name'], individuals[3]['name'])
#os.system(command)


############## Finally Running VAAST !!!


print "Start VAAST Analysis"
simulations = []

#Simulation One - Dominant, locus heterogeneity yes
simulation = {}
simulation['name'] = 'simulation_one'
simulation['parameters'] = '-iht d -lh y -fast_gp -r 0.001 --gp 10000 -m lrt  -k --codon_bias'
simulation['target'] = 'cdr/ramires.cdr'
simulations.append(simulation)

#Simulation Two - Dominant, locus heterogeneity no
simulation = {}
simulation['name'] = 'simulation_two'
simulation['parameters'] = '-iht d -lh n -fast_gp --gp 10000 -r 0.001  -m lrt  -k'
simulation['target'] = 'cdr/ramires.cdr'
simulations.append(simulation)

#Simulation Three - Dominant, penetrance complete, locus heterogenity no
simulation = {}
simulation['name'] = 'simulation_three'
simulation['parameters'] = '-iht d --penetrance c -lh n --gp 10000 -fast_gp -r 0.00035  -m lrt  -k'
simulation['target'] = 'cdr/ramires.cdr'
simulations.append(simulation)

#Simulation Four - Trio Analysis, dominant
simulation = {}
simulation['name'] = 'simulation_four'
simulation['parameters'] = '-iht d --mode lrt --codon_bias --gp 10000 -k --trio cdr/parents.cdr'
simulation['target'] = 'cdr/ramires.cdr'
simulations.append(simulation)


#Simulation Five - Trio Analysis, recessive
simulation = {}
simulation['name'] = 'simulation_five'
simulation['parameters'] = '-iht r --mode lrt --codon_bias --gp 10000 -k --trio cdr/parents.cdr'
simulation['target'] = 'cdr/ramires.cdr'
simulations.append(simulation)

#Simulation Six - -iht r  -lh n  -fast_gp  -d 1e4  -o test  -r 0.00035  -m lrt  -k
simulation = {}
simulation['name'] = 'simulation_six'
simulation['parameters'] = '-pnt c --mode lrt -iht d -lh n -fast_gp --gp 10000 --significance 2.4e-6 -r 0.00035 -m lrt -k'
simulation['target'] = 'cdr/affecteds.cdr'
simulations.append(simulation)

#Simulation Seven - Recessive, locus heterogenity no
simulation = {}
simulation['name'] = 'simulation_seven'
simulation['parameters'] = '--mode lrt -iht r -lh n -fast_gp --gp 10000 --significance 2.4e-6 -r 0.00035 -k'
simulation['target'] = 'cdr/affecteds.cdr'
simulations.append(simulation)

#Simulation 8 - 
simulation = {}
simulation['name'] = 'simulation8'
simulation['parameters'] = '-pnt c -iht d -lh n -fast_gp --gp 10000 --significance 2.4e-6 -r 0.00035 -m lrt -k'
simulation['target'] = 'cdr/affecteds.cdr'
simulations.append(simulation)



for simulation in simulations:
    command = "time %s/VAAST -mp1 8 %s -o simulations/%s %s %s %s && /lgc/programs/VAAST_Code_1.0.3/bin/vaast_tools/simple_output.pl simulations/%s.vaast > simulations/%s.genelist" % (vaast_dir, simulation['parameters'], simulation['name'], genes, background, simulation['target'], simulation['name'], simulation['name'])
    print command
    #os.system(command)
die()


#../bin/VAAST  -iht r  -lh n  -fast_gp  -d 1e4  -o test  -r 0.00035  -m lrt  -k data/easy-hg18-chr16.gff3 data/189genomes-chr16.cdr data/miller_output.cdr
#../bin/VAAST  -iht r  -lh n  -fast_gp  -d 1e6 -o test  -m lrt  -k data/easy-hg18-chr16.gff3 data/189genomes-chr16.cdr data/miller_output.cdr
#../bin/VAAST  -iht r  -lh n  -fast_gp  -d 1e4  -o test  -r 0.00035  -m lrt  -k
 

#die()

##command = "time %s/VAAST -iht r  -lh n  -fast_gp  -d 1e6 -o  %s -m lrt  -k %s %s %s.cdr" % (vaast_dir, affected_filename, genes, background, affected_filename)
#command = "time %s/VAAST -mp 4 -iht r -lh y -pnt c  -fast_gp  -m lrt -d 1e6 -k --outfile %s %s %s Exome.cdr" % (vaast_dir, affected_filename, genes, background)
##print command
#os.system(command)

#command = "~/projects/trio_analysis/bin/VAAST_Code_1.0.1/bin/vaast_tools/simple_output.pl vcf_withchr.vaast > vaast.genelist"
#os.system(command)

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


#../bin/vaast_tools/simple_output.pl 




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

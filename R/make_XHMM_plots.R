##########################################################
# First, you should make sure that the XHMM R source code
# is installed by following:
# http://atgu.mgh.harvard.edu/xhmm/download.shtml#xhmm_R_scripts
##########################################################
library(xhmmScripts)


##########################################################
# CHANGE TO DESIRED OUTPUT PATH FOR PLOTS:
##########################################################
PLOT_PATH = "plots/"


##########################################################
# CHANGE TO FULL PATH AND PREFIX FOR ALL XHMM PIPELINE OUTPUT:
#
# (e.g., DATA.RD.txt, DATA.RD_PCA.PC.txt, 
# DATA.PCA_normalized.filtered.sample_zscores.RD.txt, etc.)
#
##########################################################
JOB_PREFICES = "./DATA"


##########################################################
# CHANGE TO FULL PATH OF PLINK/Seq FORMAT EXON-TO-GENE LIST:
#
# THIS CAN BE GENERATED FOR EXOME.interval_list USING THE FOLLOWING PLINK/Seq COMMAND:
#
# pseq . loc-intersect --group refseq --locdb /path/to/locdb --file ./EXOME.interval_list --out ./annotated_targets.refseq
#
##########################################################
JOB_TARGETS_TO_GENES = "./annotated_targets.refseq.loci"



##########################################################
# CAN BE (OPTIONALLY) POPULATED FROM EITHER CHOICE BELOW:
##########################################################
SAMPLE_FEATURES = NULL


##########################################################
# CHOICE A: CHANGE TO FULL PATH OF PLINK/Seq PEDIGREE FILE
# (SEE  http://atgu.mgh.harvard.edu/plinkseq/input.shtml#ped):
##########################################################
#PEDIGREE_FILE = "/path/to/pseq_pedigree_file"
#PEDIGREE_DATA = readPedigreeFile(PEDIGREE_FILE)
#SAMPLE_FEATURES = pedigreeDataToBinarySampleProperties(PEDIGREE_DATA)


##########################################################
# CHOICE B: OR, INSTEAD OF A PLINK/Seq PEDIGREE FILE, USE A PLINK/Seq SAMPLE PHENOTYPES FILE
# (SEE  https://atgu.mgh.harvard.edu/plinkseq/input.shtml#phe):
##########################################################
#PHENOTYPES_FILE = "/path/to/pseq_phenotypes_file"
#PHENOTYPES_DATA = readPhenotypesFile(PHENOTYPES_FILE)
#SAMPLE_FEATURES = phenotypeDataToBinarySampleProperties(PHENOTYPES_DATA)



XHMM_plots(PLOT_PATH, JOB_PREFICES, JOB_TARGETS_TO_GENES, SAMPLE_FEATURES)

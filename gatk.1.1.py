#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime



__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2011, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#run example
#python gatk.py -i alignment/exome.sorted.bam

parser = OptionParser()
usage = "usage: %prog [options] -f reads1.fastq -q reads2.fastq"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted in BAM format", metavar="BAM")
                  
(options, args) = parser.parse_args()

input_file=options.input_file



#PROGRAMS
bwa_dir="/lgc/programs/bwa-0.5.9"
gatk_dir="/lgc/programs/GenomeAnalysisTK-1.1-23-g8072bd9"
pic_dir="/lgc/programs/picard-tools-1.52"
st_dir="/lgc/programs/samtools-0.1.17"


sting_dir = "/lgc/programs/Sting/dist/"

sting_res = "/lgc/programs/Sting/public/R/"

#FOLDERS

log_dir = "logs"

#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"
reference="/lgc/datasets/hg19/Homo_sapiens_assembly19.fasta"
dbsnp="/lgc/datasets/dbsnp/dbsnp-135.vcf"

omni="/lgc/datasets/gatk_data/b37_resources/1000G_omni2.5.b37.sites.vcf"
hapmap="/lgc/datasets/gatk_data/b37_resources/hapmap_3.3.b37.sites.vcf"
indels="/lgc/datasets/gatk_data/b37_resources/1000G_biallelic.indels.b37.vcf"
exome_dataset="/lgc/datasets/exome_datasets/292_illumina_samples_from_bi_and_washu.vcf"

class Gatk():
    def __init__(self, input_file):
      print "Starting GATK..."
      
      #You can change the order of the tasks :D
      starttime = datetime.datetime.now()
      print "Start Time: %s" % (starttime)
      input_file = self.AddorReplaceGroups(input_file, True)
      
      input_file = self.MarkDuplicates(input_file, True)
      
      input_file = self.LocalRealignmentAroundIndels(input_file, True)
      
      input_file = self.BaseQualityScoreRecalibration(input_file, True)
      
      input_file = self.UnifierGenotyper(input_file, True)
      die()
      self.VariantQualityScoreRecalibration(input_file)
      #self.VariantEval(input_file)
      
      timetaken = datetime.datetime.now() - starttime
      print "Time Taken: %s" % (timetaken)
      
      
      #don't forget to clean files after processing!!!!!
      
    def AddorReplaceGroups(self, input_file, flag=True):
      
      
      print "AddorReplaceGroups..."
      
      filename = self.getfilename(input_file)
      output_file = "%s.rg.bam" % (filename)
      
      command = """java -jar %s/AddOrReplaceReadGroups.jar \
      I=%s \
      O=%s \
      SORT_ORDER=coordinate \
      RGID=EXOME RGLB=EXOME RGPL=SOLID RGPU=EXOME RGSM=EXOME CREATE_INDEX=True \
      VALIDATION_STRINGENCY=LENIENT \
      """ % (pic_dir, input_file, output_file)
      
      
      if flag:
            os.system(command)

      return (output_file)
      
    def MarkDuplicates(self, input_file, flag=True):
      
	
      print "Mark Duplicates..."
      
      filename = self.getfilename(input_file)
      output_file = "%s.dedup.bam" % (filename)
      
      
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/MarkDuplicates.jar \
      INPUT=%s \
      REMOVE_DUPLICATES=true \
      VALIDATION_STRINGENCY=LENIENT \
      AS=true \
      METRICS_FILE=%s.dedup.metrics \
      OUTPUT=%s """ % (pic_dir, input_file, filename, output_file)
      
      if flag:
            os.system(command)
      
      #Reindex BAM File
      command = "%s/samtools index %s %s.bai" % (st_dir, output_file, output_file)
      if flag:
	     os.system(command)
      
      return (output_file)
      
    def LocalRealignmentAroundIndels(self, input_file, flag=True):
      
      print "Starting LocalRealignmentAroundIndels..."
      
      
      filename = self.getfilename(input_file)
      output_file = "%s.real.fixed.bam" % (filename)
      
      
      print "Starting RealignerTargetCreator..."
      command = """
      java  -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T RealignerTargetCreator \
      -l INFO \
      -R %s \
      -I %s \
      -B:dbsnp,vcf %s \
      -B:indels,vcf %s \
      -log %s/RealignerTargetCreator.log \
      -o %s.intervals """ % (gatk_dir, reference, input_file, dbsnp, indels, log_dir, filename)

      #-B:intervals,BED $EXON_CAPTURE_FILE \
      #-L $EXON_CAPTURE_FILE \
      
      if flag:
	os.system(command)
      
      print "Starting IndelRealigner..."
      # 1.2 IndelRealigner 
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T IndelRealigner \
      -R %s \
      -I %s \
      -targetIntervals %s.intervals \
      -log %s/IndelRealigner.log \
      -o %s.real.bam """ % (gatk_dir, reference, input_file, filename, log_dir, filename)
      if flag:
	os.system(command)
      
      print "FixMate"
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/FixMateInformation.jar \
      INPUT=%s.real.bam \
      OUTPUT=%s \
      SO=coordinate \
      VALIDATION_STRINGENCY=LENIENT \
      CREATE_INDEX=true
      """ % (pic_dir, filename, output_file)
      if flag:
            os.system(command)
      
      return(output_file)
    
    def BaseQualityScoreRecalibration(self, input_file, flag=True):
      print "Starting Base Quality Score Recalibration..."
      
      filename = self.getfilename(input_file)
      output_file = "%s.recal.bam" % (filename)
      
      #3.1 CountCovariates Before
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T CountCovariates \
      -l INFO \
      -I %s \
      -R %s \
      --default_read_group EXOME --default_platform SOLID \
      --solid_nocall_strategy PURGE_READ --solid_recal_mode SET_Q_ZERO_BASE_N \
      -B:dbsnp,VCF %s \
      -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \
      -recalFile %s.before.csv \
      -log %s/CountCovariates_before.log \
      -nt 4
      """ % (gatk_dir, input_file, reference, dbsnp, filename, log_dir)
      if flag:
	     os.system(command)
      # java -Xmx2g -jar GenomeAnalysisTK.jar -T BaseRecalibrator -I sample.bam -R human_g1k_v37.fasta -knownSites dbsnp_135.b37.excluding_sites_after_129.vcf -knownSites hapmap_3.3.b37.sites.vcf -knownSites 1000G_omni2.5.b37.sites.vcf -o sample.bam.recal_data.grp --covariate QualityScoreCovariate --covariate ReadGroupCovariate --covariate BinaryTagCovariate --covariate ContextCovariate --covariate CycleCovariate --solid_nocall_strategy PURGE_READ --solid_recal_mode SET_Q_ZERO_BASE_N
      


      #3.2 TableRecalibration (try -baq RECALCULATE)
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T TableRecalibration \
      -l INFO \
      -I %s \
      -R %s \
      -recalFile %s.before.csv \
      -log %s/TableRecalibration.log \
      --default_read_group EXOME_RMS --default_platform SOLID \
      --solid_nocall_strategy PURGE_READ --solid_recal_mode SET_Q_ZERO_BASE_N \
      -o %s """ % (gatk_dir, input_file, reference, filename, log_dir, output_file)
      if flag:
	     os.system(command)
      
      #Reindex Bam FILE
      command = "%s/samtools index %s %s.bai" % (st_dir, output_file, output_file)
      if flag:
	     os.system(command)
      
      # 3.2.1 CountCovariates After
      command = """
      java -Xmx12g -jar %s/GenomeAnalysisTK.jar -T CountCovariates \
      -l INFO \
      -I %s \
      -R %s \
      -B:dbsnp,VCF %s \
      --default_read_group EXOME_RMS --default_platform SOLID \
      --solid_nocall_strategy PURGE_READ --solid_recal_mode SET_Q_ZERO_BASE_N \
      -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \
      -recalFile %s.after.csv \
      -log %s/CountCovariates_after.log \
      -nt 4 """ % (gatk_dir, output_file, reference, dbsnp, filename, log_dir)
      if flag:
	     os.system(command)
      
      
      
      #3.3.2 Analyze Covariates - Generate Graphs Before  and After Table Recalibration

      #Create output folders
      if flag:
	     os.system("mkdir %s.recal.stats.before" % (filename))
      if flag:
	     os.system("mkdir %s.recal.stats.after" % (filename))
      
      # #Generate graphs before
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/AnalyzeCovariates.jar \
      -l INFO  \
      --recal_file %s.before.csv \
      -outputDir %s.recal.stats.before/ \
      -Rscript /usr/bin/Rscript \
      -ignoreQ 5 \
      -resources %s """ % (sting_dir, filename, filename, sting_res)
      if flag:
	     os.system(command)
      
       
      #Generate graphs after
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/AnalyzeCovariates.jar \
      -l INFO  \
      --recal_file %s.after.csv \
      -outputDir %s.recal.stats.after/ \
      -Rscript /usr/bin/Rscript \
      -ignoreQ 5 \
      -resources %s """ % (sting_dir, filename, filename, sting_res)
      if flag:
	     os.system(command)
      
      
      return(output_file)
    
    def UnifierGenotyper(self, input_file, flag=True):
      

      print "Starting UnifierGenotyper..."
      
      filename = self.getfilename(input_file)
      output_file = "%s.raw.vcf" % (filename)
      # # #Standard Raw VCF
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T UnifiedGenotyper \
      -l INFO \
      -I %s \
      -R %s \
      -B:dbsnp,VCF %s \
      -glm BOTH \
      -stand_call_conf 50.0 \
      -stand_emit_conf 10.0 \
      -dcov 800 \
      -A AlleleBalance \
      -A DepthOfCoverage \
      -A FisherStrand \
      -o  %s\
      -log %s/UnifiedGenotyper.log \
      -nt 4
      """ % (gatk_dir, input_file, reference, dbsnp, output_file, log_dir)
      if flag:
	     os.system(command)
      #-B:targetIntervals,BED $EXON_CAPTURE_FILE \
      #-B:intervals,BED $EXON_CAPTURE_FILE \
      return (output_file)
      
    def VariantQualityScoreRecalibration(self, input_file, flag=True):
      print "Starting Variant Quality Score Recalibration..."
      filename = self.getfilename(input_file)
      #output_file = "%s.snps.vcf" % (filename)
      
      print "SelectVariants - Select SNPS from the unified genotyper raw VCF"
      
      command = """      
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T SelectVariants \
      -R %s \
      -B:variant,VCF %s \
      -snps \
      -log %s/SelectVariants.snps.log \
      -o %s.snps.vcf """ % (gatk_dir, reference, input_file, log_dir, filename)
      if flag:
	os.system(command)
      
      #create folder
      command = "mkdir VariantRecalibrator"
      if flag:
	os.system(command)
      
      print "Start VariantRecalibrator"
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T VariantRecalibrator \
      -R %s \
      -B:input,VCF %s.snps.vcf \
      -B:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 %s \
      -B:omni,VCF,known=false,training=true,truth=false,prior=12.0 %s \
      -B:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 %s \
      -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
      --maxGaussians 4 \
      -recalFile VariantRecalibrator/%s.recal \
      -tranchesFile VariantRecalibrator/%s.tranches \
      -rscriptFile VariantRecalibrator/%s.plots.R \
      -log %s/VariantRecalibrator.log \
      -nt 4 """ % (gatk_dir, reference, filename, hapmap, omni, dbsnp, filename, filename, filename, log_dir)
      if flag:
	os.system(command)
      
      #--percentBadVariants 0.05 \
      #-mode SNP \


      #Apply Recalibrator
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T ApplyRecalibration \
      -R %s \
      -B:input,VCF %s.snps.vcf \
      -tranchesFile VariantRecalibrator/%s.tranches \
      -recalFile VariantRecalibrator/%s.recal \
      --ts_filter_level 99.0 \
      -log %s/ApplyRecalibration.log \
      -o %s.snps.recal.vcf """ % (gatk_dir, reference, filename, filename, filename, log_dir, filename)
      if flag:
	os.system(command)
      
      #Filtering SNPS Variant Filtration
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T VariantFiltration \
      -R %s \
      -B:variant,VCF %s.snps.recal.vcf \
      --clusterWindowSize 10 \
      --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
      --filterName "HARD_TO_VALIDATE" \
      --filterExpression "DP < 5 " \
      --filterName "LowCoverage" \
      --filterExpression "QUAL < 30.0 " \
      --filterName "VeryLowQual" \
      --filterExpression "QUAL > 30.0 && QUAL < 50.0 " \
      --filterName "LowQual" \
      --filterExpression "QD < 1.5 " \
      --filterName "LowQD" \
      --filterExpression "SB > -10.0 " \
      --filterName "StrandBias" \
      -o %s.snps.recal.filtered.vcf """ % (gatk_dir, reference, filename, filename)
      if flag:
	os.system(command)
      
      print "SelectVariants.Select Indels from the Unified Genotyper raw VCF"
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T SelectVariants \
      -R %s \
      -B:variant,VCF %s \
      -indels \
      -log %s/SelectVariants.indels.log \
      -o %s.indels.vcf """ % (gatk_dir, reference, input_file, log_dir, filename)
      if flag:
	os.system(command)
      
      #Variant Filtration
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T VariantFiltration \
      -R %s \
      -B:variant,VCF %s.indels.vcf \
      --filterExpression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0" \
      --filterName GATK_standard \
      --missingValuesInExpressionsShouldEvaluateAsFailing \
      -log %s/VariantFiltration.indels.log \
      -o %s.indels.filtered.vcf """ % (gatk_dir, reference, filename, log_dir, filename)
      if flag:
	os.system(command)
      
      #CombineVariants: Combine Recalibrated SNPs and Filtered  Indels
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T CombineVariants \
      -R %s \
      -B:variant,VCF %s.snps.recal.vcf \
      -B:variant,VCF %s.indels.filtered.vcf \
      -log %s/CombineVariants.log \
      -o %s.snps.indels.vcf """ % (gatk_dir, reference, filename, filename, log_dir, filename)
      if flag:
	os.system(command)
	
      command = """grep "^#" %s.snps.indels.vcf > exome.final.vcf""" % (filename)
      if flag:
	os.system(command)
	
      command = """grep "PASS" %s.snps.indels.vcf >> exome.final.vcf""" % (filename)
      if flag:
	os.system(command)
      
    def VariantEval(self, input_file, flag=True):
    
      filename = self.getfilename(input_file)
      #VariantEval
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T VariantEval \
      -R %s \
      -B:dbsnp,VCF %s \
      -o %s.eval.gatkreport \
      -B:eval,VCF $INPUT \
      -l INFO 
      -nt 4 """ % (gatk_dir, reference, dbsnp, filename, input_file)
      
      if flag:
	os.system(command)

      
    def getfilename(self, filepath):
      
      filename = ".".join(filepath.split("/")[-1].split(".")[:-1])
      return filename
      
if __name__ == '__main__': 
  Gatk(input_file)


# ##########################################STEP 2 MarkDuplicates
# 
#MARKDUPLICATES
#java -Xmx4g -jar $PIC_DIR/MarkDuplicates.jar \
#INPUT=$OUT_DIR/exome.real.bam \
#REMOVE_DUPLICATES=true \
#VALIDATION_STRINGENCY=LENIENT \
#AS=true \
#METRICS_FILE=alignment/picard/exome.dedup.metrics \
#OUTPUT=$OUT_DIR/exome.real.dedup.bam
#AS = ASSUME_SORTED
# 


# 
# #Ignore Next 2 Steps !!!!
# #2.1 AddorReplaceGroups
# java -jar $PIC_DIR/AddOrReplaceReadGroups.jar \
# I=$OUT_DIR/exome.real.dedup.bam \
# O=$OUT_DIR/exome.real.dedup.add.bam \
# SORT_ORDER=coordinate \
# RGID=EXOME RGLB=EXOME RGPL=illumina RGPU=EXOME RGSM=EXOME CREATE_INDEX=True \
# VALIDATION_STRINGENCY=LENIENT \
# TMP_DIR=/var/tmp/bio712
# # 
# # #2.2 FixMate --- should be used after realigment ???
# # java -jar  $PIC_DIR/FixMateInformation.jar \
# # INPUT=$OUT_DIR/exome.real.dedup.bam \
# # OUTPUT=$OUT_DIR/exome.real.dedup.matefixed.bam \
# # SORT_ORDER=coordinate \
# # VALIDATION_STRINGENCY=LENIENT \
# # TMP_DIR=/var/tmp/bio712
# 
# 
# ##########################################STEP 3 Base quality score recalibration
# 
# # #3.1 CountCovariates Before
# java -Xmx12g -jar $GATK_DIR/GenomeAnalysisTK.jar -T CountCovariates \
# -l INFO \
# -I $OUT_DIR/exome.real.dedup.add.bam \
# -R $REFERENCE \
# -B:dbsnp,VCF $DBSNP \
# -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \
# -recalFile $OUT_DIR/exome.recal_data.before.csv \
# -log $LOG_DIR/CountCovariates_before.log \
# -nt 4
# # 
# # #3.2 TableRecalibration (try -baq RECALCULATE)
# java -jar $GATK_DIR/GenomeAnalysisTK.jar -T TableRecalibration \
# -l INFO \
# -I $OUT_DIR/exome.real.dedup.add.bam \
# -R $REFERENCE \
# -recalFile $OUT_DIR/exome.recal_data.before.csv \
# -log $LOG_DIR/TableRecalibration.log \
# -o $OUT_DIR/exome.real.dedup.recal.bam
# 
# #Reindex Bam FILE
# $ST_DIR/samtools index $OUT_DIR/exome.real.dedup.recal.bam $OUT_DIR/exome.real.dedup.recal.bam.bai
# 
# # 3.2.1 CountCovariates After
# java -Xmx12g -jar $GATK_DIR/GenomeAnalysisTK.jar -T CountCovariates \
# -l INFO \
# -I $OUT_DIR/exome.real.dedup.recal.bam \
# -R $REFERENCE \
# -B:dbsnp,VCF $DBSNP \
# -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \
# -recalFile $OUT_DIR/exome.recal_data.after.csv \
# -log $LOG_DIR/CountCovariates_after.log \
# -nt 4
# 
# 
# ##########################################STEP 4 UnifierGenotyper
# # 
# # #Standard Raw VCF
# java -Xmx15g -jar $GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper \
# -l INFO \
# -I $OUT_DIR/exome.real.dedup.recal.bam \
# -R $REFERENCE \
# -B:intervals,BED $EXON_CAPTURE_FILE \
# -B:dbsnp,VCF $DBSNP \
# -glm BOTH \
# -stand_call_conf 50.0 \
# -stand_emit_conf 20.0 \
# -dcov 300 \
# -A AlleleBalance \
# -A DepthOfCoverage \
# -A FisherStrand \
# -o $OUT_DIR/exome.raw.vcf \
# -log $LOG_DIR/UnifiedGenotyper.log \
# -nt 4
# #-B:targetIntervals,BED $EXON_CAPTURE_FILE \
# 
# ##########################################STEP 5 Variant quality score recalibration
# 
# #SelectVariants.Select SNPS from the unified genotyper raw VCF
# java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar -T SelectVariants \
# -R $REFERENCE \
# -B:variant,VCF $OUT_DIR/exome.raw.vcf \
# -snps \
# -log $LOG_DIR/SelectVariants.snps.log \
# -o $OUT_DIR/exome.snps.raw.vcf
# 
# 
# # #create folder
# mkdir $OUT_DIR/VariantRecalibrator
# #VariantRecalibrator
# java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar -T VariantRecalibrator \
# -R $REFERENCE \
# -B:input,VCF $OUT_DIR/exome.snps.raw.vcf \
# -B:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $HAPMAP \
# -B:omni,VCF,known=false,training=true,truth=false,prior=12.0 $OMNI \
# -B:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $DBSNP \
# -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
# --maxGaussians 4 \
# -recalFile $OUT_DIR/VariantRecalibrator/exome.recal \
# -tranchesFile $OUT_DIR/VariantRecalibrator/exome.tranches \
# -rscriptFile $OUT_DIR/VariantRecalibrator/exome.plots.R \
# -log $LOG_DIR/VariantRecalibrator.log \
# -nt 4

#--percentBadVariants 0.05 \
#-mode SNP \

# 
#  
# #Apply Recalibrator
# java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar -T ApplyRecalibration \
# -R $REFERENCE \
# -B:input,VCF $OUT_DIR/exome.snps.raw.vcf \
# -tranchesFile $OUT_DIR/VariantRecalibrator/exome.tranches \
# -recalFile $OUT_DIR/VariantRecalibrator/exome.recal \
# --ts_filter_level 99.0 \
# -log $LOG_DIR/ApplyRecalibration.log \
# -o $OUT_DIR/exome.recal.snps.vcf

#Filtering SNPS Variant Filtration
# java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar -T VariantFiltration \
# -R $REFERENCE \
# -B:variant,VCF $OUT_DIR/exome.recal.snps.vcf \
# --clusterWindowSize 10 \
# --filterExpression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
# --filterName "HARD_TO_VALIDATE" \
# --filterExpression "DP < 5 " \
# --filterName "LowCoverage" \
# --filterExpression "QUAL < 30.0 " \
# --filterName "VeryLowQual" \
# --filterExpression "QUAL > 30.0 && QUAL < 50.0 " \
# --filterName "LowQual" \
# --filterExpression "QD < 1.5 " \
# --filterName "LowQD" \
# --filterExpression "SB > -10.0 " \
# --filterName "StrandBias" \
# -o $OUT_DIR/exome.snps.recal.filtered.vcf



# #SelectVariants.Select Indels from the Unified Genotyper raw VCF
# java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar -T SelectVariants \
# -R $REFERENCE \
# -B:variant,VCF $OUT_DIR/exome.raw.vcf \
# -indels \
# -log $LOG_DIR/SelectVariants.indels.log \
# -o $OUT_DIR/exome.indels.raw.vcf
## 
## #Variant Filtration
#java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar -T VariantFiltration \
#-R $REFERENCE \
#-B:variant,VCF $OUT_DIR/exome.indels.raw.vcf \
#--filterExpression "QD < 2.0 || ReadPosRankSum < -20.0 || FS > 200.0" \
#--filterName GATK_standard \
#--missingValuesInExpressionsShouldEvaluateAsFailing \
#-log $LOG_DIR/VariantFiltration.indels.log \
#-o $OUT_DIR/exome.indels.filtered.vcf
## 
## #CombineVariants: Combine Recalibrated SNPs and Filtered  Indels
#java -Xmx4g -jar $GATK_DIR/GenomeAnalysisTK.jar -T CombineVariants \
#-R $REFERENCE \
#-B:variant,VCF $OUT_DIR/exome.recal.snps.vcf \
#-B:variant,VCF $OUT_DIR/exome.indels.filtered.vcf \
#-log $LOG_DIR/CombineVariants.log \
#-o $OUT_DIR/exome.snps.indels.vcf

#Include Deph Of Coverage!!

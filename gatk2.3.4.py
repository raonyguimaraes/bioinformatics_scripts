#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime



__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2013, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "2.3.4"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#run example
#python gatk.py -i alignment/exome.sorted.bam

parser = OptionParser()
usage = "usage: %prog [options] -i reads.bam"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted in BAM format", metavar="BAM")
parser.add_option("-t", dest="target_array",
                  help="Target Array", metavar="BEDFILE")
                  
(options, args) = parser.parse_args()

input_file=options.input_file
target_array=options.target_array


#PROGRAMS

gatk_dir="/lgc/programs/gatk/dist"
gatk_dir='/lgc/programs/GenomeAnalysisTK-2.3-4-g57ea19f'
gatk_dir='/lgc/programs/GenomeAnalysisTK-2.3-5-g49ed93c'


pic_dir="/lgc/programs/picard-tools-1.82"
st_dir="/lgc/programs/samtools"


sting_dir = "/lgc/programs/gatk/dist"
sting_res = "/lgc/programs/gatk/public/R/"
#/lgc/programs/gatk/public/R/scripts/org/broadinstitute/sting/gatk/walkers/variantrecalibration/plot_Tranches.R
#FOLDERS

log_dir = "logs"
reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"
#reference="/lgc/datasets/hg19/Homo_sapiens_assembly19.fasta"
dbsnp="/lgc/datasets/dbsnp/dbsnp-135.vcf"
dbsnp="/lgc/datasets/dbsnp/137/00-All.vcf"

omni="/lgc/datasets/gatk_data/b37_resources/1000G_omni2.5.b37.sites.vcf"
hapmap="/lgc/datasets/gatk_data/b37_resources/hapmap_3.3.b37.sites.vcf"
indels="/lgc/datasets/gatk_data/b37_resources/1000G_phase1.indels.b37.vcf"
exome_dataset="/lgc/datasets/exome_datasets/292_illumina_samples_from_bi_and_washu.vcf"

exon_10bp = "/lgc/datasets/ucsc/refseq_plus10.withoutchr.bed"


class Gatk():
    def __init__(self, input_file):
      print "Starting GATK..."
      
      #You can change the order of the tasks :D
      starttime = datetime.datetime.now()
      print "Start Time: %s" % (starttime)
      
      input_file = self.MarkDuplicates(input_file, False)

      input_file = self.LocalRealignmentAroundIndels(input_file, False)
      
      input_file = self.BaseQualityScoreRecalibration(input_file, True)
      die()

      #CALCULATE Depth of Coverage
      self.DepthofCoverage(input_file, True)
      input_file = self.UnifierGenotyper(input_file, True)
      self.VariantQualityScoreRecalibration(input_file)
      #die()
      #self.VariantEval(input_file)
      #self.Dindel(input_file)
      
      timetaken = datetime.datetime.now() - starttime
      print "Time Taken: %s" % (timetaken)
      
      
      #don't forget to clean files after processing!!!!!
      
      
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
      CREATE_INDEX=true \
      METRICS_FILE=%s.dedup.metrics \
      OUTPUT=%s """ % (pic_dir, input_file, filename, output_file)
      if flag:
	     os.system(command)
      
      
      
      
      
      #Reindex BAM File
      command = "%s/samtools index %s %s.bai" % (st_dir, output_file, output_file)
      # if flag:
            # os.system(command)

      #perl -e 'print "@RG\tID:ga\tSM:hs\tLB:ga\tPL:Illumina\n@RG\tID:454\tSM:hs\tLB:454\tPL:454\n"' > rg.txt
      #samtools merge -rh rg.txt - ga.bam 454.bam | samtools rmdup - - | samtools rmdup -s - aln.bam
      # bwa samse -r @RG\tID:IDa\tSM:SM\tPL:Illumina ref.fa my.sai my.fastq > my.sam
      #java -jar AddOrReplaceReadGroups I=my.bam O=myGr.bam LB=whatever PL=illumina PU=whatever SM=whatever


      return (output_file)
      
    def LocalRealignmentAroundIndels(self, input_file, flag=True):
      
      print "Add or Replace Read Groups..."

      # command = 'java -jar AddOrReplaceReadGroups I=my.bam O=myGr.bam LB=whatever PL=illumina PU=whatever SM=whatever' 
      
      filename = self.getfilename(input_file)
      
      output_file = "%s.rg.bam" % (filename)

      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/AddOrReplaceReadGroups.jar \
      I=%s \
      O=%s \
      LB=exome \
      PL=solid \
      VALIDATION_STRINGENCY=LENIENT \
      PU=exome \
      SM=exome \
      """ % (pic_dir, input_file, output_file)
      
      

      if flag:
            os.system(command)
           

      #Reindex BAM File
      print 'Reindexing'
      command = "%s/samtools index %s %s.bai" % (st_dir, output_file, output_file)
      if flag:
            os.system(command)

      print "Starting LocalRealignmentAroundIndels..."
      
      
      filename = self.getfilename(output_file)
      output_file = "%s.real.fixed.bam" % (filename)
      
      
      print "Starting RealignerTargetCreator..."
      command = """
      java  -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T RealignerTargetCreator \
      -l INFO \
      -R %s \
      -I %s.bam \
      --known %s \
      --known %s \
      -nt 8 \
      --fix_misencoded_quality_scores \
      -log %s/RealignerTargetCreator.log \
      -o %s.intervals """ % (gatk_dir, reference, filename, dbsnp, indels, log_dir, filename)

      #-B:intervals,BED $EXON_CAPTURE_FILE \
      #-L $EXON_CAPTURE_FILE \
      
      if flag:
            os.system(command)
      # die()

      
      print "Starting IndelRealigner..."
      # 1.2 IndelRealigner 
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T IndelRealigner \
      -R %s \
      -I %s.bam \
      --fix_misencoded_quality_scores \
      -targetIntervals %s.intervals \
      -log %s/IndelRealigner.log \
      -o %s.real.bam """ % (gatk_dir, reference, filename, filename, log_dir, filename)
      
      
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
      
      # #3.1 CountCovariates Before
      # command = """
      # java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T CountCovariates \
      # -l INFO \
      # -I %s \
      # -R %s \
      # -knownSites %s \
      # -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \
      # -recalFile %s.before.csv \
      # -log %s/CountCovariates_before.log \
      # -nt 4
      # """ % (gatk_dir, input_file, reference, dbsnp, filename, log_dir)
      
      #3.1 BaseRecalibrator Before
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T BaseRecalibrator \
      -l INFO \
      -I %s \
      -R %s \
      -knownSites %s \
      -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov ContextCovariate \
      --filter_mismatching_base_and_quals \
      -o recal_data.grp \
      """ % (gatk_dir, input_file, reference, dbsnp)
      #-log %s/BaseRecalibrator.log \
      #, log_dir
      if flag:
            os.system(command)
      
      #
      #3.2 Print Reads
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T PrintReads \
      -l INFO \
      -I %s \
      -R %s \
      -BQSR recal_data.grp \
      -log %s/TableRecalibration.log \
      -o %s """ % (gatk_dir, input_file, reference, log_dir, output_file)
      if flag:
            os.system(command)

      
      
      #Reindex Bam FILE
      command = "%s/samtools index %s %s.bai" % (st_dir, output_file, output_file)
      if flag:
	     os.system(command)
      
 #      # 3.2.1 CountCovariates After
 #      command = """
 #      java -Xmx12g -jar %s/GenomeAnalysisTK.jar -T CountCovariates \
 #      -l INFO \
 #      -I %s \
 #      -R %s \
 #      -knownSites %s \
 #      -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate \
 #      -recalFile %s.after.csv \
 #      -log %s/CountCovariates_after.log \
 #      -nt 4 """ % (gatk_dir, output_file, reference, dbsnp, filename, log_dir)
 #      if flag:
	# os.system(command)
      
      
      
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
    
    def DepthofCoverage(self, input_file, flag=True):
      print "Starting DepthOfCoverage..."
      filename = self.getfilename(input_file)
      
      command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T DepthOfCoverage \
      -I %s \
      -R %s \
      -o %s_depthofcoverage \
      -L %s \
      -log %s/DepthofCoverage.log \
      """ % (gatk_dir, input_file, reference, filename, target_array, log_dir)
      if flag:
	os.system(command)
	
      #[-geneList refSeq.sorted.txt] \
      #-ct 4 -ct 6 -ct 10 \
      #[-L my_capture_genes.interval_list]
    
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
      --dbsnp %s \
      --genotype_likelihoods_model BOTH \
      -stand_call_conf 50.0 \
      -stand_emit_conf 10.0 \
      -dcov 50 \
      -A AlleleBalance \
      -A DepthOfCoverage \
      -A FisherStrand \
      -o  %s\
      -log %s/UnifiedGenotyper.log \
      -L %s \
      -nt 4
      """ % (gatk_dir, input_file, reference, dbsnp, output_file, log_dir, target_array)
      if flag:
	os.system(command)
      #-B:targetIntervals,BED $EXON_CAPTURE_FILE \
      #-B:intervals,BED $EXON_CAPTURE_FILE \
      #-L targets.interval_list
      return (output_file)
      
    def VariantQualityScoreRecalibration(self, input_file, flag=True):
      print "Starting Variant Quality Score Recalibration..."
      filename = self.getfilename(input_file)
      #output_file = "%s.snps.vcf" % (filename)
      
      print "SelectVariants - Select SNPS from the unified genotyper raw VCF"
      
      command = """      
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T SelectVariants \
      -R %s \
      --variant %s \
      -selectType SNP \
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
      -input %s.snps.vcf \
      -resource:hapmap,known=false,training=true,truth=true,prior=15.0 %s \
      -resource:omni,known=false,training=true,truth=false,prior=12.0 %s \
      -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 %s \
      -resource:exome_dataset,known=true,training=true,truth=true,prior=8.0 %s \
      -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ \
      --maxGaussians 4 \
      -recalFile VariantRecalibrator/%s.recal \
      -tranchesFile VariantRecalibrator/%s.tranches \
      -rscriptFile VariantRecalibrator/%s.plots.R \
      -log %s/VariantRecalibrator.log \
      -nt 4 """ % (gatk_dir, reference, filename, hapmap, omni, dbsnp, exome_dataset, filename, filename, filename, log_dir)
      if flag:
	os.system(command)
      
      #--percentBadVariants 0.05 \
      #-mode SNP \


      #Apply Recalibrator
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T ApplyRecalibration \
      -R %s \
      -input %s.snps.vcf \
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
      --variant %s.snps.recal.vcf \
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
      --variant %s \
      -selectType INDEL \
      -log %s/SelectVariants.indels.log \
      -o %s.indels.vcf """ % (gatk_dir, reference, input_file, log_dir, filename)
      if flag:
	os.system(command)
      
      #Variant Filtration
      command = """
      java -Xmx4g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T VariantFiltration \
      -R %s \
      --variant %s.indels.vcf \
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
      --variant %s.snps.recal.vcf \
      --variant %s.indels.filtered.vcf \
      -genotypeMergeOptions UNIQUIFY \
      -log %s/CombineVariants.log \
      -o %s.snps.indels.vcf """ % (gatk_dir, reference, filename, filename, log_dir, filename)
      if flag:
	os.system(command)
	
      command = """grep "^#" %s.snps.indels.vcf > exome.passed.vcf""" % (filename)
      if flag:
	os.system(command)
	
      command = """grep "PASS" %s.snps.indels.vcf >> exome.passed.vcf""" % (filename)
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

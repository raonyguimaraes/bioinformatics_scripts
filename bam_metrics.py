#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
		
#run example
#python bam_metrics.py -i exome.bam

parser = OptionParser()


parser.add_option("-i", dest="input_file",
                  help="BAM File in BAM format", metavar="BAM")
                  
(options, args) = parser.parse_args()

input_file=options.input_file


#variables
pic_dir="/lgc/programs/picard-tools-1.52"
reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"

exon = "../../input/bed/Design_Annotation_files/Target_Regions/SeqCap_EZ_Exome_v2.bed"

out_dir="alignment_metrics"
os.system("mkdir alignment_metrics")


class Bam_metrics():
    def __init__(self, input_file):
      print "Starting Bam_metrics..."
      #piccard tools
      filename = ".".join(input_file.split("/")[-1].split(".")[:-1])
      
      
      #CollectAlignmentSummaryMetrics
      command = """
      java -jar -Xmx40g %s/CollectAlignmentSummaryMetrics.jar \
      I=%s \
      O=%s.AlignmentSummaryMetrics \
      R=%s \
      VALIDATION_STRINGENCY=LENIENT """ % (pic_dir, input_file, filename, reference)
      os.system(command)
      
      #CollectGcBiasMetrics
      command = """
      java -jar %s/CollectGcBiasMetrics.jar \
      R=%s \
      I=%s \
      O=%s.b37_1kg.GcBiasMetrics \
      CHART=%s.b37_1kg.GcBiasMetrics.pdf \
      VALIDATION_STRINGENCY=LENIENT """ % (pic_dir, reference, input_file, filename, filename)
      os.system(command)

      #CollectInsertSizeMetrics
      command = """
      java -jar %s/CollectInsertSizeMetrics.jar \
      I=%s \
      O=%s.b37_1kg.CollectInsertSizeMetrics \
      H=%s.b37_1kg.CollectInsertSizeMetrics.pdf \
      VALIDATION_STRINGENCY=LENIENT """ % (pic_dir, input_file, filename, filename)
      os.system(command)
      
      #MeanQualityByCycle
      command = """
      java -jar %s/MeanQualityByCycle.jar \
      I=%s \
      O=%s.b37_1kg.MeanQualityByCycle \
      CHART=%s.b37_1kg.MeanQualityByCycle.pdf \
      VALIDATION_STRINGENCY=LENIENT """ % (pic_dir, input_file, filename, filename)
      os.system(command)
      
      #QualityScoreDistribution
      command = """
      java -jar %s/QualityScoreDistribution.jar \
      I=%s \
      O=%s.b37_1kg.QualityScoreDistribution \
      CHART=%s.b37_1kg.QualityScoreDistribution.pdf \
      VALIDATION_STRINGENCY=LENIENT """ % (pic_dir, input_file, filename, filename)
      os.system(command)
      
      #BamIndexStats
      command = """
      java -jar %s/BamIndexStats.jar \
      INPUT=%s \
      VALIDATION_STRINGENCY=LENIENT """ % (pic_dir, input_file)
      os.system(command)
      
      
      
      
      #CalculateHsMetrics WholeGenome or CalculateHsMetrics ????
      command = """
      java -jar -Xmx4g %s/CalculateHsMetrics.jar \
      INPUT=%s \
      OUTPUT=%s.b37_1kg.HsMetrics \
      BAIT_INTERVALS=%s \
      TARGET_INTERVALS=%s \
      VALIDATION_STRINGENCY=LENIENT """ % (pic_dir, input_file, filename, exon, exon)
      #os.system(command)
      
      
      
      #Samstat
      
      


if __name__ == '__main__': 
  Bam_metrics(input_file)

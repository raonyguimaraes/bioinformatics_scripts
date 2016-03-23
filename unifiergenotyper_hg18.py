#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
from time import time
import datetime
import re
#Example
#python /projects/12exomes/data/unifiergenotyper.py -i /projects/12exomes/data/Exome_3_EDS.realigned-recalibrated.bam /projects/12exomes/data/Exome_4_ELS.realigned-recalibrated.bam /projects/12exomes/data/Exome_5_LS.realigned-recalibrated.bam /projects/12exomes/data/Exome_6_DC.realigned-recalibrated.bam -o quartet_only_SNPs_.vcf -t /lgc/datasets/exome_targets/SureSelect_All_Exon_V2_hg19.20110105.bed


__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2012, The Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "1.0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Production"
		
#run example
#python unifiergenotyper.py -i /projects/12exomes/data/Exome_3_EDS.realigned-recalibrated.bam /projects/12exomes/data/Exome_4_ELS.realigned-recalibrated.bam /projects/12exomes/data/Exome_5_LS.realigned-recalibrated.bam /projects/12exomes/data/Exome_6_DC.realigned-recalibrated.bam

parser = OptionParser()
usage = "usage: %prog [options] -f reads1.fastq -q reads2.fastq"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_file",
                  help="BAM File Sorted in BAM format", metavar="BAM")#, nargs=12
parser.add_option("-t", dest="target_array",
                  help="Target Array", metavar="BEDFILE")
parser.add_option("-o", dest="output",
                  help="Output Filename", metavar="VCFFILE")
                  
(options, args) = parser.parse_args()

input_file=options.input_file



filename = ".".join(input_file.split("/")[-1].split(".")[:-1])

target_array=options.target_array

#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta"
#reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.fasta"
#dbsnp="/lgc/datasets/dbsnp/dbsnp-135chr.vcf"

reference="/lgc/datasets/gatk_data/b37/human_g1k_v37chr.fasta"
reference="/lgc/datasets/hg18/all/all.fasta"
#dbsnp="/lgc/datasets/dbsnp/137/00-All.vcf"

# gatk_dir="/lgc/programs/GenomeAnalysisTK-2.3-9-ge5ebf34"
# gatk_dir="/lgc/programs/GenomeAnalysisTK-1.1-23-g8072bd9"
gatk_dir="/lgc/programs/GenomeAnalysisTK-2.5-2-gf57256b"
pic_dir="/lgc/programs/picard-tools-1.82"


#--dbsnp %s \

#add replace groups
command = """java -jar %s/AddOrReplaceReadGroups.jar \
      I=%s \
      O=%s.rg.bam \
      SORT_ORDER=coordinate \
      RGID=EXOME RGLB=EXOME RGPL=ILLUMINA RGPU=EXOME RGSM=EXOME CREATE_INDEX=True \
      VALIDATION_STRINGENCY=LENIENT \
      """ % (pic_dir, input_file, filename)
#os.system(command)

command = """
      java -Xmx40g -Djava.io.tmpdir=/projects/tmp -jar %s/GenomeAnalysisTK.jar -T UnifiedGenotyper \
      -l INFO \
      -I %s \
      -R %s \
      -o  %s\
      -L %s \
      --downsampling_type NONE \
      -nt 4
      """ % (gatk_dir, input_file, reference, options.output, re.escape(target_array))
#--fix_misencoded_quality_scores \
#      -rf ReassignMappingQuality -DMQ 255 -dcov 50 \
      
#--fix_misencoded_quality_scores \
#-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60
      #-rf ReassignMappingQuality -DMQ 255 \
      
      #--output_mode EMIT_ALL_CONFIDENT_SITES \
#grep "241371863" child.less.vcf 
#java -Xmx12g -jar GenomeAnalysisTK.jar -nt 6 -T UnifiedGenotyper -I $1 -o $1.SNP.vcf -R ref.fasta -glm SNP -metrics $1.SNP.metrics -rf MappingQuality -mmq 10 -rf BadCigar java -Xmx8g -jar GenomeAnalysisTK.jar -R re#f.fasta -T CombineVariants --variant *.SNP.vcf --variant *.INDEL.vcf -o output.CombinedVariants.vcf

print command
os.system(command)    


#vcf_file = open("%s.vcf" % (filename), 'r')
#vcf_output = open(filename, 'w')
#for line in vcf_file :
    #if line.startswith("#"):
	#vcf_output.write(line)
    #else:
	#vcf_output.write(line.replace('chr', ''))
#vcf_output.close()
#vcf_file.close()

#command = 'mv %s %s.vcf' % (filename, filename)
#os.system(command)


      
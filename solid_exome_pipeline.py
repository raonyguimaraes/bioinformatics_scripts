#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os

__author__ = "Raony Guimarãees"
__copyright__ = "Copyright 2012, Solid Exome Pipeline"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Development"
		
#run example
#python bwa.py -f ../../input/exome/Ot729_index7_Exome_Sergio.Pena_RDNP_07262011-reads1-110716_I123_FCC043NABXX_L2_index7_1.fastq -q ../../input/exome/Ot729_index7_Exome_Sergio.Pena_RDNP_07262011-reads2-110716_I123_FCC043NABXX_L2_index7_2.fastq

parser = OptionParser()
usage = "usage: %prog [options] -f reads1.fastq -q reads2.fastq"
parser = OptionParser(usage=usage)

parser.add_option("-f", dest="reads1",
                  help="reads 1 in FASTQ format", metavar="FASTQ")
parser.add_option("-r", dest="reads2",
                  help="reads 2 in FASTQ format", metavar="FASTQ")
#parser.add_option("-t", dest="target_array",
                  #help="Target Array", metavar="BEDFILE")
                  
(options, args) = parser.parse_args()


reads1=options.reads1
reads1_name = reads1.split('/')[-1]
reads2=options.reads2
reads2_name = reads2.split('/')[-1]


#target_array=options.target_array

#print "Alignment"
#command = 'python /lgc/scripts/bwa.py -f %s -q %s' % (reads1, reads2)
#os.system(command)
#print "GATK"
#command = 'python /lgc/scripts/gatk.py -i exome.sorted.bam -t %s' % (target_array)
#os.system(command)

#BWA Index Normal
#command = "/lgc/programs/bwa-0.5.9/bwa index -a bwtsw /lgc/datasets/gatk_data/b37/normal/human_g1k_v37.fasta"
#BWA Index ColorSpace
#command = "/lgc/programs/bwa-0.5.9/bwa index -c -a bwtsw /lgc/datasets/gatk_data/b37/colorspace/human_g1k_v37.fasta"


bwa_dir = "/lgc/programs/bwa-0.5.9"
bfast_dir = "/lgc/programs/bfast_bwa/bfast/bfast"
human_reference = "/lgc/datasets/gatk_data/b37/bfast/human_g1k_v37.fasta"

#r3 = "/lgc/datasets/exomes/hpf/ngs_analyses/work/llau/brazil/samples/Exome_10_EF-20120504095627-gatk1.1.28f-hg19/read-convert-R3-20120504095627-gatk1.1.28f/Exome_10_EF.1.fastq"
#f3 = "/lgc/datasets/exomes/hpf/ngs_analyses/work/llau/brazil/samples/Exome_10_EF-20120504095627-gatk1.1.28f-hg19/read-convert-F3-20120504095627-gatk1.1.28f/Exome_10_EF.1.fastq"

#father_r3_dir = "/lgc/datasets/exomes/hpf/ngs_analyses/work/llau/brazil/samples/Exome_11_AP-20120504101053-gatk1.1.28f-hg19/read-convert-R3-20120504101053-gatk1.1.28f/"
#father_f3_dir = "/lgc/datasets/exomes/hpf/ngs_analyses/work/llau/brazil/samples/Exome_11_AP-20120504101053-gatk1.1.28f-hg19/read-convert-F3-20120504101053-gatk1.1.28f/"

#dirList=os.listdir(father_r3_dir)
#for fname in dirList:
    #print fname


#F3 = 75bp
#R3 = 35bp

#Aligning with bwa for Reverse
#command = "%s/bwa aln -c -t 8  %s %s > Exome_10_EF.1.R.sai 2> Exome_10_EF.1.R.log" % (bwa_dir, human_reference, r3)
#print command
#os.system(command)

#command = "%s/bwa samse %s Exome_10_EF.1.R.sai %s > out.sam" % (bwa_dir, human_reference, r3)
#os.system(command)


#BFAST for FORWARD F3
#CREATE A REFERENCE GENOME
def convert_reference():
  command = "%s/bfast fasta2brg -f %s " % (bfast_dir, human_reference)
  os.system(command)
  command = "%s/bfast fasta2brg -f %s -A 1" % (bfast_dir, human_reference)
  os.system(command)
  

def index_bfast():
  #CREATE INDEX FOR GENOME REFERENCE
  print 'CREATE INDEX FOR GENOME REFERENCE'
  #50base pairs
  masks = ['1111111111111111111111', '111110100111110011111111111', '10111111011001100011111000111111', '1111111100101111000001100011111011', '111111110001111110011111111', '11111011010011000011000110011111111', '1111111111110011101111111', '111011000011111111001111011111', '1110110001011010011100101111101111', '111111001000110001011100110001100011111']
  #25 base pairs
  alternative_masks = ['111111111111111111', '11110100110111101010101111', '11111111111111001111', '1111011101100101001111111', '11110111000101010000010101110111', '1011001101011110100110010010111', '1110110010100001000101100111001111', '1111011111111111111', '11011111100010110111101101', '111010001110001110100011011111']
  index = 1
  for mask in masks:
    print 'index %s' % (index)
    command = "time %s/bfast index -n 8 -f %s -m %s -w 14 -i %s -A 1 2>output.%s.log" % (bfast_dir, human_reference, mask, index, index)
    print command
    #os.system(command)
    index += 1
    
  #index for bwa
  command = "%s/bfast bwtindex -a bwtsw -A 1 %s" % (bfast_dir, human_reference)
  os.system(command)

def alignment():
  #For 75bp F3 Forward
  #MATCH
  command = "time %s/bfast match -n 16 -f %s -A 1 -r %s > %s.bmf 2>match.%s.log" % (bfast_dir, human_reference, reads1, reads1_name, reads1_name)
  os.system(command)
  
  #die()
  #LOCALALIGN
  #command = "%s/bfast localalign -n 8 -f %s -m bfast.matches.fasta.reads.bmf -A 1 > bfast.aligned.fasta.reads.baf" % (bfast_dir, human_reference)
  #os.system(command)
  #command = "%s/bfast postprocess -n 8 -f %s -i bfast.aligned.fasta.reads.baf -A 1 > bfast.reported.fasta.reads.sam" % (bfast_dir, human_reference)
  #os.system(command)

  #For 35bp R3 Reverse
  command = "time %s/bfast bwaaln -t 16 -c %s %s > %s.bmf 2>bwaalign.%s.log" % (bfast_dir, human_reference, reads2, reads2_name, reads2_name)
  os.system(command)
  
  #bring them together:
  command = "time %s/bfast localalign -f %s -1 %s.bmf -2 %s.bmf -A 1 -t -U -n 16 > %s.%s.baf 2>localalign.%s.%s.log" % (bfast_dir, human_reference, reads1_name, reads2_name, reads1_name, reads2_name, reads1_name, reads2_name)
  os.system(command)

  #But the postprocess step, which was done in a few minutes for single end, can take > 100 hours on 16 CPUs for 50 Mio read pairs:
  command = "time %s/bfast postprocess -f %s -i %s.%s.baf -a 3 -A 1 -z -t -n 16 >%s.%s.sam" % (bfast_dir, human_reference, reads1_name, reads2_name, reads1_name, reads2_name)
  #print command
  os.system(command)
  #samtools view -bS reads.sam > reads.bam   
  #samtools sort reads.bam aln-sorted
  #samtools index aln-sorted.bam
  
  #command = "time %s/bfast postprocess -f %s -i reads.baf -a 3 -A 1 -R -z -v 160 -s 20 -S 4.0 > reads.sam 2>postprocess.log" % (bfast_dir, human_reference)
   
  
  
print 'reference'
#convert_reference() # Done
print 'indexing'
#index_bfast()
alignment()

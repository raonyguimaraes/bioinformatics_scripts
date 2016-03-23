#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
import csv

parser = OptionParser()
parser.add_option("-v", dest="vcf_file",
                  help="VCF File to Annotate", metavar="VCF")
#parser.add_option("-o", dest="vcf_output",
                  #help="VCF File to Annotate", metavar="VCF")

(options, args) = parser.parse_args()



#"Create a folder and get inside..."
##os.chdir( path ). 
#print os.chdir('.')
#remove CHR from file
command = "python /lgc/scripts/vcf_without_chr.py -i %s" % (options.vcf_file)
os.system(command)


#Annotate with SNPEFF
filename = os.path.splitext(options.vcf_file)[0]

#Sort vcf file before annotation
os.system("grep '^#' %s.vcf > output.vcf" % (filename))
os.system("grep -E -v '^X|^Y|^MT|^#|^GL' %s.vcf | sort -n -k1 -k2 >> output.vcf" % (filename))
os.system("grep -E '^X' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))
os.system("grep -E '^Y' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))
os.system("grep -E '^MT' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))
#Ignore mutations from non-reference
#os.system("grep -E '^GL' %s.vcf | sort -k1,1d -k2,2n >> output.vcf" % (filename))

os.system("mv output.vcf %s.vcf" % (filename))


#Create SnpEff annotations
command = "java -Xmx4G -jar /lgc/programs/snpeff_2.0.5/snpEff_2_0_5/snpEff.jar eff -c /lgc/programs/snpeff_2.0.5/snpEff_2_0_5/snpEff.config -v -onlyCoding true -i vcf -o vcf GRCh37.64 %s > snpEff_output.vcf" % (options.vcf_file)

#command = "java -Xmx4g -jar /lgc/programs/snpEff_3_0/snpEff.jar eff -c /lgc/programs/snpEff_3_0/snpEff.config -v GRCh37.66 -o gatk %s > snpEff_output.vcf" % (options.vcf_file)
os.system(command)

#removing snpeff versions
#command = " sed -i 's/SnpEffVersion=\"2.0.5d/SnpEffVersion=\"2.0.5 /g' snpEff_output.vcf"
#os.system(command)
#command = " sed -i 's/SnpEffVersion=\"SnpEff 3.0f /SnpEffVersion=\"2.0.5 /g' snpEff_output.vcf"
#os.system(command)


#all in one big command!
command = "java -Xmx10G -jar /lgc/programs/GenomeAnalysisTK-1.6-13-g91f02df/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R /lgc/datasets/gatk_data/b37/human_g1k_v37.fasta \
-A SnpEff \
--alwaysAppendDbsnpId \
--variant %s \
--snpEffFile snpEff_output.vcf \
-L %s \
-o %s.snpeff.dbsnp.exomeserver.vcf \
--dbsnp /lgc/datasets/dbsnp/137/00-All.vcf \
--resource:dbsnp137 /lgc/datasets/dbsnp/137/00-All.vcf \
--resource:exome_server /lgc/datasets/exome_variation_server/ESP6500.vcf \
--resource:1000genomes /lgc/datasets/1000genomes/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf \
-E 1000genomes.AF \
-E exome_server.AAC \
-E exome_server.GS \
-E exome_server.PP \
-E exome_server.CDP \
-E exome_server.MAF \
-E exome_server.CA \
-E exome_server.TAC \
-E dbsnp137.dbSNPBuildID \
-E dbsnp137.GMAF \
-E dbsnp137.PM \
" % (options.vcf_file, options.vcf_file, filename)

#command = "java -Xmx10G -jar /lgc/programs/GenomeAnalysisTK-2.0-35-g2d70733/GenomeAnalysisTK.jar \
#-T VariantAnnotator \
#-R /lgc/datasets/gatk_data/b37/human_g1k_v37.fasta \
#-A SnpEff \
#--alwaysAppendDbsnpId \
#--variant %s \
#--snpEffFile snpEff_output.vcf \
#-L %s \
#-o %s.snpeff.dbsnp.exomeserver.vcf \
#--dbsnp /lgc/datasets/dbsnp/137/00-All.vcf \
#--resource:dbsnp137 /lgc/datasets/dbsnp/137/00-All.vcf \
#--resource:exome_server /lgc/datasets/exome_variation_server/ESP6500.vcf \
#--resource:1000genomes /lgc/datasets/1000genomes/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf \
#-E 1000genomes.AF \
#-E exome_server.AAC \
#-E exome_server.GS \
#-E exome_server.PP \
#-E exome_server.CDP \
#-E exome_server.MAF \
#-E exome_server.CA \
#-E exome_server.TAC \
#-E dbsnp137.dbSNPBuildID \
#-E dbsnp137.GMAF \
#-E dbsnp137.PM \
#" % (options.vcf_file, options.vcf_file, filename)

os.system(command)

##removing SNPEFF using VCFTOOLS
os.environ["PERL5LIB"] = "/lgc/programs/vcftools/lib/perl5/site_perl/"
command =  "cat %s.snpeff.dbsnp.exomeserver.vcf | /lgc/programs/vcftools/bin/vcf-annotate -r INFO/SNPEFF_AMINO_ACID_CHANGE,INFO/SNPEFF_CODON_CHANGE,INFO/SNPEFF_GENE_BIOTYPE,INFO/SNPEFF_TRANSCRIPT_ID,INFO/SNPEFF_EXON_ID > %s.snpeff.dbsnp.exomeserver.cutted.vcf" % (filename, filename)
os.system(command)

snpeff_tags = ['##INFO=<ID=SNPEFF_AMINO--snpEffFile snpEff_output.vcf -L exome.sorted.dedup.real.fixed.recal.raw.vcf -o exome.sorted.dedup.real.fixed.recal.raw.snpeff.dbsnp.exomeserver.vcf --dbsnp /lgc/datasets/dbsnp/137/00-All.vcf --resource:dbsnp /lgc/datasets/dbsnp/137/00-All.vcf --resource:exome_server /lgc/datasets/exome_variation_server/ESP6500.vcf --resource:1000genomes /lgc/datasets/1000genomes/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf -E 1000genomes.AF -E dbsnp.GMAF -E dbsnp.CLN -E dbsnp.dbSNPBuildID -E exome_server.AAC -E exome_server.GS -E exome_server.PP -E exome_server.CDP -E exome_server.MAF -E exome_server.CA -E exome_server.TAC_ACID_CHANGE', '##INFO=<ID=SNPEFF_CODON_CHANGE', '##INFO=<ID=SNPEFF_GENE_BIOTYPE', '##INFO=<ID=SNPEFF_TRANSCRIPT_ID', '##INFO=<ID=SNPEFF_EXON_ID']

#Annotate with annovar!!!!
#command = "/lgc/scripts/annovar.py -i %s.vcf"  % (filename)
#os.system(command)
ann_dir="/lgc/programs/annovar"

os.system('rm -rf annovar')
os.system('mkdir annovar')

command = "%s/convert2annovar.pl --format vcf4 --includeinfo %s.vcf > annovar/%s.annovar" % (ann_dir, filename, filename)
os.system(command)

#--ver1000g 1000g2011may
command = "%s/summarize_annovar.pl --ver1000g 1000g2012apr -verdbsnp 135 --buildver hg19 annovar/%s.annovar %s/humandb -outfile annovar/%s" % (ann_dir, filename, ann_dir, filename)
os.system(command)


##annotate with SIFT WEB
#command = 'python /lgc/scripts/sift.py -i %s' % (options.vcf_file)
#os.system(command)

##annotate with Annovar
#command = 'python /lgc/scripts/polyphen2.py -i %s' % (options.vcf_file)
#os.system(command)


##index annotated variants
##semaphonre for both for saving time
##read files into index to include into vcf
#sift_variants = {}
#pp2_variants = {}



annovar_file_path =  "annovar/%s.genome_summary.csv" % (filename)
annovar_file = csv.reader(open(annovar_file_path, "rb"))
header = annovar_file.next()
#print header

variants = {}
for line in annovar_file:
    variant = {}
    sift_index = header.index('AVSIFT')
    polyphen_index = header.index('LJB_PolyPhen2')
    #Extra Tags from Annovar
    variant['conserved'] = line[header.index('Conserved')].replace(';', '|')
    variant['segdup'] = line[header.index('SegDup')]
    variant['esp6500']= line[header.index('ESP6500_ALL')]
    variant['avsift'] = line[header.index('AVSIFT')]
    variant['ljb_phylop'] = line[header.index('LJB_PhyloP')]
    variant['ljb_phylop_pred'] = line[header.index('LJB_PhyloP_Pred')]
    variant['ljb_sift'] = line[header.index('LJB_SIFT')]
    variant['ljb_sift_pred'] = line[header.index('LJB_SIFT_Pred')]
    variant['ljb_polyphen2'] = line[header.index('LJB_PolyPhen2')]
    variant['ljb_polyphen2_pred'] = line[header.index('LJB_PolyPhen2_Pred')]
    variant['ljb_lrt'] = line[header.index('LJB_LRT')]
    variant['ljb_lrt_pred'] = line[header.index('LJB_LRT_Pred')]
    variant['ljb_mutationtaster'] = line[header.index('LJB_MutationTaster')]
    variant['ljb_mutationtaster_pred'] = line[header.index('LJB_MutationTaster_Pred')]
    variant['ljb_gerp'] = line[header.index('LJB_GERP++')]

    
    variant['sift'] = line[sift_index]
    variant['polyphen'] = line[polyphen_index]
    
    chr_index = header.index('Otherinfo')
    pos_index = chr_index+1
    
    variant_id = "%s-%s" % (line[chr_index], line[pos_index])
    variants[variant_id] = variant 
   
#print len(variants)
vcf_file=open('%s.snpeff.dbsnp.exomeserver.vcf' % (filename), 'r')
#out_file=open('%s' % (options.vcf_output), 'w')
vcf_full=open('%s.fullannotation.vcf' % (filename), 'w')

for line in vcf_file:
	if line.startswith('#'):
	    info_tag = line.split(',')[0]
	    if info_tag not in snpeff_tags:
		#out_file.writelines(line)
		vcf_full.writelines(line)
	else:
	    line = line.split('\t')
	    variant_id = "%s-%s" % (line[0], line[1])
	    if variant_id in variants:
		if variants[variant_id]['sift'] != '':
		    line[7] = line[7]+';ANN_SIFT=%s' % (variants[variant_id]['sift'])
		if variants[variant_id]['polyphen'] != '':
		    line[7] = line[7]+';ANN_POLYPHEN=%s' % (variants[variant_id]['polyphen'])
		#extra tags
		string_full_annotation = line[7]
		if variants[variant_id]['conserved'] != '':
		    string_full_annotation = string_full_annotation+';ANN_CONSERVED=%s' % (variants[variant_id]['conserved'])
		if variants[variant_id]['segdup'] != '':
		    string_full_annotation = string_full_annotation+';ANN_SEGDUP=%s' % (variants[variant_id]['segdup'])
		if variants[variant_id]['esp6500'] != '':
		    string_full_annotation = string_full_annotation+';ANN_ESP6500=%s' % (variants[variant_id]['esp6500'])
		if variants[variant_id]['avsift'] != '':
		    string_full_annotation = string_full_annotation+';ANN_AVSIFT=%s' % (variants[variant_id]['avsift'])
		if variants[variant_id]['ljb_phylop'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_PHYLOP=%s' % (variants[variant_id]['ljb_phylop'])
		if variants[variant_id]['ljb_phylop_pred'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_PHYLOP_PRED=%s' % (variants[variant_id]['ljb_phylop_pred'])
		if variants[variant_id]['ljb_sift'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_SIFT=%s' % (variants[variant_id]['ljb_sift'])
		if variants[variant_id]['ljb_sift_pred'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_SIFT_PRED=%s' % (variants[variant_id]['ljb_sift_pred'])
		if variants[variant_id]['ljb_polyphen2'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_POLYPHEN2=%s' % (variants[variant_id]['ljb_polyphen2'])
		if variants[variant_id]['ljb_polyphen2_pred'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_POLYPHEN2_PRED=%s' % (variants[variant_id]['ljb_polyphen2_pred'])
		if variants[variant_id]['ljb_lrt'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_LRT=%s' % (variants[variant_id]['ljb_lrt'])
		if variants[variant_id]['ljb_lrt_pred'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_LRT_PRED=%s' % (variants[variant_id]['ljb_lrt_pred'])
		if variants[variant_id]['ljb_mutationtaster'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_MUTATIONTASTER=%s' % (variants[variant_id]['ljb_mutationtaster'])
		if variants[variant_id]['ljb_mutationtaster_pred'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_MUTATIONTASTER_PRED=%s' % (variants[variant_id]['ljb_mutationtaster_pred'])
		if variants[variant_id]['ljb_gerp'] != '':
		    string_full_annotation = string_full_annotation+';ANN_LJB_GERP=%s' % (variants[variant_id]['ljb_gerp'])
		
		    
		#out_file.writelines("\t".join(line))
		#add all annotations to line
		line[7] = string_full_annotation
		vcf_full.writelines("\t".join(line))
		
	    else:
		#out_file.writelines("\t".join(line))
		vcf_full.writelines("\t".join(line))
#out_file.close()
vcf_full.close()


###cleaning files
#command = "rm -f snpEff*" 
#os.system(command)
#command = "rm -f *.snpeff.vcf*"
#os.system(command)
#os.system("rm -f *.idx")
#command = "rm -f %s.snpeff.dbsnp.exomeserver.vcf" % (filename)
#os.system(command)
#command = "rm -f %s.snpeff.dbsnp.exomeserver.cutted.vcf" % (filename)
##do not remove to fix the bug
#os.system(command)
##remove temp files from annovar
#command = 'mv annovar/%s.genome_summary.csv .' % (filename, ) 
#os.system(command)
#os.system('rm -rf annovar')

#command = "rm -rf annovar" % (filename)
#os.system(command)

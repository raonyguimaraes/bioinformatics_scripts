from optparse import OptionParser
import os


parser = OptionParser()
parser.add_option("-v", dest="vcf_file",
                  help="VCF File to Annotate using dbSNP135", metavar="VCF")

(options, args) = parser.parse_args()

#Annotate with SNPEFF

filename = options.vcf_file.split('/')[-1].split('.')[0]
print filename
#die()
#Create SnpEff annotations
command = "java -Xmx4G -jar /lgc/programs/snpeff/snpEff_2_0_5/snpEff.jar eff -c /lgc/programs/snpeff/snpEff_2_0_5/snpEff.config -v -onlyCoding true -i vcf -o vcf GRCh37.64 %s > snpEff_output.vcf" % (options.vcf_file)
os.system(command)

#all in one big command!
command = "java -Xmx4G -jar ~/lgc/programs/GenomeAnalysisTK-1.4-30-gf2ef8d1/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R ~/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta \
-A SnpEff \
--variant %s \
--snpEffFile snpEff_output.vcf \
-L %s \
-o %s.snpeff.dbsnp135.exomeserver.vcf \
--dbsnp /home/raony/lgc/datasets/dbsnp/dbsnp-135.vcf \
--resource:dbsnp135 /home/raony/lgc/datasets/dbsnp/dbsnp-135.vcf \
--resource:exome_server /home/raony/lgc/datasets/dbsnp/dbsnp-135.vcf \
-E dbsnp135.GENEINFO \
-E dbsnp135.GMAF \
-E dbsnp135.SAO \
-E dbsnp135.SSR \
-E dbsnp135.SCS \
-E dbsnp135.WGT \
-E dbsnp135.VC \
-E dbsnp135.CLN \
-E dbsnp135.PM \
-E dbsnp135.TPA \
-E dbsnp135.PMC \
-E dbsnp135.S3D \
-E dbsnp135.SLO \
-E dbsnp135.NSF \
-E dbsnp135.NSM \
-E dbsnp135.HD \
-E dbsnp135.OM \
-E dbsnp135.G5A \
-E dbsnp135.G5 \
-E dbsnp135.dbSNPBuildID \
-E exome_server.TAC \
-E exome_server.AA_AC \
-E exome_server.EA_AC \
-E exome_server.PH \
-E exome_server.FG \
-E exome_server.AAC \
-E exome_server.PP" % (options.vcf_file, options.vcf_file, filename)

os.system(command)

#annotate a vcf with dbsnp135 ids
#command = "java -Xmx4G -jar ~/lgc/programs/GenomeAnalysisTK-1.4-30-gf2ef8d1/GenomeAnalysisTK.jar \
#-T VariantAnnotator \
#-R ~/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta \
#--variant %s \
#-L %s \
#-o %s.dbsnp135.vcf \
#--dbsnp /home/raony/lgc/datasets/dbsnp/dbsnp-135.vcf" % (options.vcf_file, options.vcf_file, options.vcf_file)

#os.system(command)


#annotate vcf file with another vcf_file
#command = "java -Xmx4G -jar ~/lgc/programs/GenomeAnalysisTK-1.4-30-gf2ef8d1/GenomeAnalysisTK.jar \
#-T VariantAnnotator \
#-R ~/lgc/datasets/gatk_data/b37/human_g1k_v37.fasta \
#--variant %s \
#-L %s \
#-o %s.full_annotation.vcf \
#--resource:dbsnp135 /home/raony/lgc/datasets/dbsnp/dbsnp-135.vcf \
#--resource:exome_server /home/raony/lgc/datasets/dbsnp/dbsnp-135.vcf \
#-E dbsnp135.GENEINFO \
#-E dbsnp135.GMAF \
#-E dbsnp135.SAO \
#-E dbsnp135.SSR \
#-E dbsnp135.SCS \
#-E dbsnp135.WGT \
#-E dbsnp135.VC \
#-E dbsnp135.CLN \
#-E dbsnp135.PM \
#-E dbsnp135.TPA \
#-E dbsnp135.PMC \
#-E dbsnp135.S3D \
#-E dbsnp135.SLO \
#-E dbsnp135.NSF \
#-E dbsnp135.NSM \
#-E dbsnp135.HD \
#-E dbsnp135.OM \
#-E dbsnp135.G5A \
#-E dbsnp135.G5 \
#-E dbsnp135.dbSNPBuildID \
#-E exome_server.TAC \
#-E exome_server.AA_AC \
#-E exome_server.EA_AC \
#-E exome_server.PH \
#-E exome_server.FG \
#-E exome_server.AAC \
#-E exome_server.PP \
#-E exome_server.PP \
#" % (options.vcf_file, options.vcf_file, options.vcf_file)

annovar_file = csv.reader(open('/home/raony/projects/exome_analysis/output/combined_annovar/merge/merge.exome_summary.csv', "rb"))
header = annovar_file.next()
print header

variants = {}
for line in annovar_file:
    variant = {}
    variant['sift'] = line[10]
    variant['polyphen'] = line[11]
    variant_id = "%s-%s" % (line[15], line[16])
    variants[variant_id] = variant
    
    
print len(variants)
vcf_file=open('%s.snpeff.dbsnp135.exomeserver.vcf' % (filename), 'r')
out_file=open('%s.snpeff.dbsnp135.exomeserver.annovar.vcf' % (filename), 'w')

for line in vcf_file:
	if line.startswith('#'):
	    out_file.writelines(line)
	else:
	    line = line.split('\t')
	    variant_id = "%s-%s" % (line[0], line[1])
	    if variant_id in variants:
		if variants[variant_id]['sift'] != '':
		    line[7] = line[7]+';ANN_SIFT=%s' % (variants[variant_id]['sift'])
		if variants[variant_id]['polyphen'] != '':
		    line[7] = line[7]+';ANN_POLYPHEN=%s' % (variants[variant_id]['polyphen'])
		out_file.writelines("\t".join(line))
	    else:
		out_file.writelines("\t".join(line))
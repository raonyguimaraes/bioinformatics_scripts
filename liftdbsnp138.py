
# command = "java -jar 
import os
gatk_dir = '/lgc/programs/GenomeAnalysisTK-2.6-5-gba531bd'

command = "java -Xmx40G -jar %s/GenomeAnalysisTK.jar \
        -T VariantAnnotator \
        -R /lgc/datasets/gatk_data/b37/human_g1k_v37.fasta \
        --alwaysAppendDbsnpId \
        --variant 1000genomes.vcf.bgzip.gz \
        -L 1000genomes.vcf.bgzip.gz \
        -o 1000genomes.dbsnp138.vcf.gz \
        --dbsnp /lgc/datasets/dbsnp/138/00-All.vcf \
        --log_to_file gatkreport.dbsnp138.log \
		" % (gatk_dir)
os.system(command)

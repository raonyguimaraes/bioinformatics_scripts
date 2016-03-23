#python script to run xhmm
#sudo apt-get install liblapack-dev liblapack3

import os, subprocess

gatk_dir = '/home/raony/bin/gatk'
xhmm_dir = '/home/raony/bin/xhmm/statgen-xhmm-7ce9df09e19d'
plinkseq_dir = '/home/raony/bin/plinkseq'
# plinkseq_dir = '/projects/xhmm_analysis/bin/plinkseq-0.10'

group1 = 'ampliseq.group1.list'
group2 = 'ampliseq.group2.list'
group3 = 'ampliseq.group3.list'

exome_interval = '/home/raony/datasets/beds/AmpliSeqExome.20141113.designed_effective.bed'
human_reference = '/home/raony/datasets/b37/human_g1k_v37chr.fasta'

bam_files = [group1,group2,group3]

resources = '/home/raony/datasets/xhmm/'
import pipes

#requirement: index bam files
def index_bams():
	for group in bam_files:
		print 'group', group
		bam_file_reader = open(group, 'r')
		for bam_file in bam_file_reader:
			bam_file = pipes.quote(bam_file.strip())
			print 'bam_file', bam_file
			index_file = "%s.bai" % (bam_file)
			if not os.path.isfile(index_file):
				index_command = 'samtools index %s' % (bam_file)
				os.system(index_command)
#index_bams()
# die()
command = 'java -Xmx40g -jar %s/GenomeAnalysisTK.jar \
-T DepthOfCoverage -I %s -L %s \
-R %s \
-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
--includeRefNSites \
--countType COUNT_FRAGMENTS \
-o group1.DATA' % (gatk_dir, group1, exome_interval, human_reference)
#print command
#os.system(command)
#die()

command = 'java -Xmx40g -jar %s/GenomeAnalysisTK.jar \
-T DepthOfCoverage -I %s -L %s \
-R %s \
-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
--includeRefNSites \
--countType COUNT_FRAGMENTS \
-o group2.DATA' % (gatk_dir, group2, exome_interval, human_reference)
#print command
# os.system(command)

command = 'java -Xmx40g -jar %s/GenomeAnalysisTK.jar \
-T DepthOfCoverage -I %s -L %s \
-R %s \
-dt BY_SAMPLE -dcov 5000 -l INFO --omitDepthOutputAtEachBase --omitLocusTable \
--minBaseQuality 0 --minMappingQuality 20 --start 1 --stop 5000 --nBins 200 \
--includeRefNSites \
--countType COUNT_FRAGMENTS \
-o group3.DATA' % (gatk_dir, group3, exome_interval, human_reference)
#print command
#die()
#os.system(command)
#die()

#Combines GATK Depth-of-Coverage outputs for multiple samples (at same loci):
command = '%s/build/execs/xhmm --mergeGATKdepths -o ./DATA.RD.txt \
--GATKdepths group1.DATA.sample_interval_summary \
--GATKdepths group2.DATA.sample_interval_summary \
--GATKdepths group3.DATA.sample_interval_summary' % (xhmm_dir)
#os.system(command)
#die()


#run GATK to calculate the per-target GC content and
#create a list of the targets with extreme GC content:
command='java -Xmx40g -jar %s/GenomeAnalysisTK.jar \
-T GCContentByInterval -L %s \
-R %s \
-o ./DATA.locus_GC.txt' % (gatk_dir, exome_interval, human_reference)
#os.system(command)

#working!


command = """cat ./DATA.locus_GC.txt | awk '{if ($2 < 0.2 || $2 > 0.8) print $1}' > ./extreme_gc_targets.txt """
os.system(command)
#die()

#instal csh sudo apt-get install csh
#run PLINK/Seq to calculate the fraction of repeat-masked bases in each target and create a list of those to filter out:
command = "%s/sources/scripts/interval_list_to_pseq_reg %s > ./EXOME.targets.reg" % (xhmm_dir, exome_interval)
# print command
os.system(command)


#git clone git@bitbucket.org:statgen/plinkseq.git
#sudo apt-get install plink
command = "%s/pseq . loc-load --locdb ./EXOME.targets.LOCDB --file ./EXOME.targets.reg --group targets --out ./EXOME.targets.LOCDB.loc-load" % (plinkseq_dir)
os.system(command)


command = """%s/pseq . --noweb loc-stats --locdb ./EXOME.targets.LOCDB --group targets --seqdb %s/seqdb | \
awk '{if (NR > 1) print $_}' | sort -k1 -g | awk '{print $10}' | paste %s - | \
awk '{print $1"\t"$2}' \
> ./DATA.locus_complexity.txt """ % (plinkseq_dir, plinkseq_dir, exome_interval)
os.system(command)
# die()

command = "cat ./DATA.locus_complexity.txt | awk '{if ($2 > 0.25) print $1}' \
> ./low_complexity_targets.txt"
os.system(command)
# die()

# Filters samples and targets and then mean-centers the targets:
command = "%s/xhmm --matrix -r ./DATA.RD.txt --centerData --centerType target \
-o ./DATA.filtered_centered.RD.txt \
--outputExcludedTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
--outputExcludedSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 25 --maxMeanSampleRD 200 \
--maxSdSampleRD 150" % (xhmm_dir)

os.system(command)
#--excludeTargets ./extreme_gc_targets.txt --excludeTargets ./low_complexity_targets.txt \


# Runs PCA on mean-centered data:
command ="%s/xhmm --PCA -r ./DATA.filtered_centered.RD.txt --PCAfiles ./DATA.RD_PCA" %(xhmm_dir)
os.system(command)
# die()

# Normalizes mean-centered data using PCA information:
command = "%s/xhmm --normalize -r ./DATA.filtered_centered.RD.txt --PCAfiles ./DATA.RD_PCA \
--normalizeOutput ./DATA.PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7" % (xhmm_dir)
os.system(command)


# Filters and z-score centers (by sample) the PCA-normalized data:
command = "%s/xhmm --matrix -r ./DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
-o ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
--outputExcludedTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--outputExcludedSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
--maxSdTargetRD 30" % (xhmm_dir)
os.system(command)


# Filters original read-depth data to be the same as filtered, normalized data:
command = "%s/xhmm --matrix -r ./DATA.RD.txt \
--excludeTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
--excludeTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--excludeSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o ./DATA.same_filtered.RD.txt" % (xhmm_dir)
os.system(command)


# Discovers CNVs in normalized data:
command = "%s/xhmm --discover -p ./params.txt \
-r ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ./DATA.same_filtered.RD.txt \
-c ./DATA.xcnv -a ./DATA.aux_xcnv -s ./DATA" % (xhmm_dir)
os.system(command)


# NOTE: A description of the .xcnv file format can be found below.

# Genotypes discovered CNVs in all samples:
command = "%s/xhmm --genotype -p ./params.txt \
-r ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ./DATA.same_filtered.RD.txt \
-g ./DATA.xcnv -F %s \
-v ./DATA.vcf" % (xhmm_dir, human_reference)
os.system(command)
# die()

# And, the R visualization plots:

# Annotate exome targets with their corresponding genes:
command = "%s/pseq . loc-intersect --group refseq --locdb ./RefSeq.locdb --file ./EXOME.interval_list --out ./annotated_targets.refseq" % (plinkseq_dir)
os.system(command)

#command from plots.R
command = "%s/pseq . loc-intersect --group refseq --locdb %s/locdb --file %s --out ./annotated_targets.refseq" % (plinkseq_dir, resources, exome_interval)
os.system(command)


# Plot the XHMM pipeline and CNV discovered:
# [ First, ensure that the R scripts are installed as described here ]

# Rscript example_make_XHMM_plots.R



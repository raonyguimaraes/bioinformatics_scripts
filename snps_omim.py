# snps = open('/lgc/datasets/omim/snps_2_omim.txt', 'r')

# snp_list = []
# for line in snps:
# 	row = line.split(' ')
# 	snp_list.append(row[0])


# #now get frequencies in big file
# g1000file = open('/projects/1000genomes/integrated_call_sets/ALL.wgs.integrated_phase1_v3.20101123.snps_indels_sv.sites.vcf', 'r')


# for line in g1000file:
# 	if not line.startswith('#'):
# 		row = line.split('\t')
# 		if row[2] in snp_list:
# 			# print line
# 			info = row[7].split(';')
# 			for tag in info:
# 				if tag.startswith('AF'):
# 					print tag

		

#using vcftools

#generate a list of snps


#vcftools --vcf file1.vcf --out snps_omim.txt --recode --snps snps_rsid_omim.txt
#parse clinvar
clinvar = open('/lgc/datasets/dbsnp/138/clinvar_00-latest.vcf', 'r')
for line in clinvar:
	if not line.startswith('#'):
		# print line
		row = line.split('\t')
		info = row[7].split(';')
		for tag in info:
			if tag == 'OM':
				print line,

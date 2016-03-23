clinvar = open('/lgc/datasets/dbsnp/141/clinvar_20140702.vcf')

variants = 0
variants_before_129 = 0
pathogenic_before_129 = 0
for line in clinvar:
	if not line.startswith('#'):
		variants += 1
		row = line.split('\t')
		info = row[7].split(';')
		#get dbsnpbuild
		dbsnpbuild = 0
		for item in info:
			if item.startswith('dbSNPBuildID'):
				dbsnpbuild = int(item.split('=')[1])
				# print 'dbsnpbuild', dbsnpbuild
				if dbsnpbuild <= 129:
					variants_before_129 += 1
					#check if it's pathogenic
			if item.startswith('CLNSIG'):
				# print item
				significance = item.split('=')[1]
				# 5 is pathogenic
				# print 'significance', significance
				if dbsnpbuild <= 129 and significance == '5':
					pathogenic_before_129 += 1
					print line



print 'total variants: ', variants
print 'total variants_before_129: ', variants_before_129
print 'total pathogenic_before_129: ', pathogenic_before_129

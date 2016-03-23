fileplink = open('/projects/relatedness/18people/plink.frq', 'r')
fileplink.next()
excludefile = open('/projects/relatedness/18people/removesnps', 'w')
for line in fileplink:
	variant = line.split()
	if variant[2] == 'A' and variant[3] == 'T':
		excludefile.writelines(variant[1]+'\n')
	if variant[2] == 'T' and variant[3] == 'A':
		excludefile.writelines(variant[1]+'\n')
	if variant[2] == 'G' and variant[3] == 'C':
		excludefile.writelines(variant[1]+'\n')
	if variant[2] == 'C' and variant[3] == 'G':
		excludefile.writelines(variant[1]+'\n')

excludefile.close()

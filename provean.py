from optparse import OptionParser


#open vcf file
#submit to provean or run provean with vcf

parser = OptionParser()

parser.add_option('-i', help='vcf', dest='vcffile') #, nargs=2
(options, args) = parser.parse_args()

vcffile=options.vcffile
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])

readfile = open(vcffile, 'r')
outfile = open('provean_variants', 'w')

for line in readfile:
	if not line.startswith('#'):
		variant = line.split('\t')
		variant[0] =  variant[0].replace('chr', '')
		outputline = '%s,%s,%s,%s\n' % (variant[0], variant[1], variant[3], variant[4])
		outfile.writelines(outputline)
outfile.close()
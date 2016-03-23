
#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
import subprocess
import time


parser = OptionParser()

parser.add_option("-i", dest="vcffile",
                  help="VCF File", metavar="VCF")
(options, args) = parser.parse_args()
vcffile=options.vcffile
# control=options.control
filename = ".".join(vcffile.split("/")[-1].split(".")[:-1])
outfile = open('%s.noy.vcf' % (filename), 'w')
for line in open(vcffile, 'r'):
	if line.startswith('#'):
		outfile.writelines(line)
	else:
		row = line.split('\t')
		if not line.startswith('Y'):
			outfile.writelines("\t".join(row))

outfile.close()

# command = 'grep "^#" %s.vcf > %s.sorted.vcf' % (filename,filename)
# os.system(command)
# for i in range(1,23):
# 	command = 'grep "^%s\t" %s.vcf >> %s.sorted.vcf' % (i, filename,filename)
# 	os.system(command)


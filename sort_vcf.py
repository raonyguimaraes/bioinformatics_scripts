
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

command = 'grep "^#" %s.vcf > %s.sorted.vcf' % (filename,filename)
os.system(command)
for i in range(1,23):
	command = 'grep "^%s\t" %s.vcf >> %s.sorted.vcf' % (i, filename,filename)
	os.system(command)

command = 'grep "^X\t" %s.vcf >> %s.sorted.vcf' % (filename,filename)
os.system(command)	
command = 'grep "^Y\t" %s.vcf >> %s.sorted.vcf' % (filename,filename)
os.system(command)	
command = 'grep "^MT\t" %s.vcf >> %s.sorted.vcf' % (filename,filename)
os.system(command)	


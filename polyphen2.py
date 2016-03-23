#!/usr/bin/env python
# -*- coding: utf-8 -*-

from optparse import OptionParser
import os
import csv
import time

parser = OptionParser()
parser.add_option("-i", dest="vcf_file",
                  help="VCF File to Annotate", metavar="VCF")
                  
(options, args) = parser.parse_args()

output = open('polyphen2_variants.txt', 'w')
for line in open(options.vcf_file):
    if not line.startswith('#'):
	line = line.split('\t')
	output.writelines('chr%s:%s %s/%s\n' % (line[0], line[1], line[3], line[4]))
                  
command = "curl -F _ggi_project=PPHWeb2 -F _ggi_origin=query -F _ggi_target_pipeline=1 -F MODELNAME=HumDiv -F UCSCDB=hg19 -F SNPFUNC=m -F NOTIFYME=raonyguimaraes@gmail.com -F _ggi_batch_file=@%s -D - http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi > polyphen.response" % ('polyphen2_variants.txt')
os.system(command)
for line in open('polyphen.response'):
    if line.startswith('Set-Cookie:'):
	line = line.split(' ')
	polyid = line[1].replace('polyphenweb2=', '').replace(';', '')
#polyid = "2ae4ae351a5ef03c90e14611aa3b7ccb430c9e63"
#check results
h = os.system("wget -q http://genetics.bwh.harvard.edu/ggi/pph2/%s/1/completed.txt" % (polyid))
print 'clear'
print h

while h != 0:
    time.sleep(60)
    h = os.system("wget -q http://genetics.bwh.harvard.edu/ggi/pph2/%s/1/completed.txt" % (polyid))
    print '.',

#downloading results
os.system("wget -q http://genetics.bwh.harvard.edu/ggi/pph2/%s/1/pph2-full.txt -O pph2-full.txt" % (polyid))
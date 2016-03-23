#!/usr/bin/env python
# -*- coding: utf-8 -*-

morbidmap = open('/lgc/datasets/omim/morbidmap', 'r') 

count_mendelian_diseases = 0
for line in morbidmap:
  line_start = line
  line = line.split('|')
  genes = line[1].split(',')
  if len(genes) == 1:
    if not line[0].startswith('{'):
      if not line[0].startswith('['):
	print line_start,
	count_mendelian_diseases += 1
	
	
print 'numero de doencas mendelianas: %s' % (count_mendelian_diseases)
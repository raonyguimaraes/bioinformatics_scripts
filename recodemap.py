mapfile = open('/projects/relatedness/hgdp/HGDP_hg19.map', 'r')
newmapfile = open('/projects/relatedness/hgdp/HGDP_hg19.withchr.map', 'w')
for line in mapfile:
	row = line.strip().split('\t')
	row[1] = "chr%s:%s" % (row[0], row[-1])
	newmapfile.writelines('\t'.join(row)+'\n')

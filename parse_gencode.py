gencode = open('/lgc/datasets/gencode.v12.annotation.gtf', 'r')
for line in gencode:
  if not line.startswith('##'):
    #print line
    row = line.split('\t')
    #print row
    gene_info = row[8].split(';')[4].strip().replace('gene_name ', '').replace('"', '')
    #print gene_info
    print "\t".join([row[0].replace('chr', ''), row[3], row[4], gene_info])
    
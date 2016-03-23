#!/usr/bin/env python

#get all genes from a pathway using kegg
#http://www.biostars.org/post/show/6224/gene-pathway-association-file-for-kegg/

from SOAPpy import WSDL

wsdl = 'http://soap.genome.jp/KEGG.wsdl'
serv = WSDL.Proxy(wsdl)

#GET GENES BY PATHWAY
print 'hedgehog signaling pathway'
results = serv.get_genes_by_pathway('path:hsa04340')
print len(results)
print 'notch signaling pathway'
results = serv.get_genes_by_pathway('path:hsa04330')
print len(results)
print '	TGF-beta signaling pathway'
results = serv.get_genes_by_pathway('path:hsa04350')
print len(results)
die()
print 'Genes'
for gene_id in results:
    gene = serv.bget(gene_id)
    for line in gene.split('\n'):
	if line.startswith('NAME'):
	    print line.replace('NAME        ', '')
	    #gene_names = line.split(', ')
	    #gene_names[0] = gene_names[0].replace('NAME        ', '')
	    #print ", ".join(gene_names)
	    
#print len(results)
#print serv.bconv(' '.join(results))

#for line in result.split('\n'):
    
#result = serv.bconv('hsa:122011')
#print result


#TGF-beta signaling pathway - Homo sapiens (human) hsa04350 (PITX2)


#GET PATHWAYS BY GENE


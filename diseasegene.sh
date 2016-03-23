curl -s "ftp://anonymous:raonyguimaraes%40gmail%2Ecom@grcf.jhmi.edu/OMIM/mim2gene.txt" |\
egrep -v "#" | cut -d ' ' -f 4 | egrep -v '^\-$' |\
sort | uniq > list1.txt


mysql -N --user=genome --host=genome-mysql.cse.ucsc.edu -A  -D hg19 -e 'select  distinct
  G.geneSymbol,
  S.name
from snp132 as S,
kgXref as G,
knownGene as K where
    S.chrom=K.chrom and
    S.chromStart>=K.txStart and
    S.chromEnd<=K.txEnd and
    K.name=G.kgId 
#    /* AND something to restrict the result to YOUR list of SNPs or gene */
' | sort -t '    ' -k1,1 > list2.txt



join -1 1 -2 1 list1.txt list2.txt

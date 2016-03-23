fasta = open('human_g1k_v37.fasta')
out_fasta = open('human_g1k_v37.chr.fasta','w')

for line in fasta:
    if line.startswith('>'):
        line = line.replace('>', '>chr')
    out_fasta.writelines(line)

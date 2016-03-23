#detect which target was used for this exome

bam_files = open('/storage/received_exomes/exomes_proton_mm.list')

import pysam

outfile = open('mm_bam_files_ampliseq.list', 'w')

for bam_file in  bam_files:
	print bam_file
	bam_file = bam_file.strip()
	samfile = pysam.AlignmentFile(bam_file, "rb")
	header = samfile.header
	if 'PU' in header['RG'][0]:
		print header['RG'][0]['PU']
		if 'IonXpress' in header['RG'][0]['PU']:
			print 'achou AmpliSeq!'
			outfile.writelines(bam_file+'\n')
outfile.close()	
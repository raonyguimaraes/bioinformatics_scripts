import os


pic_dir="/storage2/programs/piccard_tools/picard-tools-1.119"
reference="/home/raony/datasets/hg19/hg19.fasta"
gatk_dir="/storage2/programs/gatk"
target = '/home/raony/datasets/beds/refseq_plus_10bp.trimmed.bed'

if not os.path.exists('rg'):
    os.makedirs('rg')

if not os.path.exists('vcfs'):
    os.makedirs('vcfs')

bam_files = open('bam_files.txt')

for bam in bam_files:
    bam = bam.strip()[2:]
    sample_id = os.path.splitext(bam)[0]
    #rg
    command =  """
        java -jar -Xmx40g -Djava.io.tmpdir=/tmp -jar %s/AddOrReplaceReadGroups.jar I=%s O=rg/%s.rg.bam LB=exome PL=illumina PU=exome SM=%s VALIDATION_STRINGENCY=SILENT
        """ % (pic_dir, bam, sample_id, sample_id)      
    print command
    os.system(command)
    os.chdir('rg')
    #index
    command = 'samtools index %s.rg.bam' % (sample_id)
    os.system(command)
    os.chdir('..')
    #call
    command = """
      java -Xmx40g -Djava.io.tmpdir=/tmp -jar %s/GenomeAnalysisTK.jar -T UnifiedGenotyper \
      -R %s \
      -I rg/%s.rg.bam \
      -l INFO \
      -stand_call_conf 50.0 \
      -stand_emit_conf 10.0 \
      -dcov 200 \
      -o  vcfs/%s.vcf \
      -log %s-UnifiedGenotyper.log \
      -L %s \
      -nt 4
      """ % (gatk_dir, reference, sample_id, sample_id, sample_id, target)
    print command
    os.system(command)
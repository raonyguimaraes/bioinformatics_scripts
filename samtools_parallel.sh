
sample=$1
describer=$(echo ${sample} | sed 's/.bam//')  
   
# Convert file from SAM to BAM format  
# samtools view -b $sample > ${describer}.uns.bam  
   
# Sort BAM file  
samtools sort $sample ${describer}.sorted
   
# index the bam file  
samtools index ${describer}.sorted.bam
   
# Revove intermediate files  
# rm ${describer}.uns.bam  
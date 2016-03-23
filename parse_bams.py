import os

bamfile = open('bam_files.txt')

paths = []

files = []
count_files = 0
for line in bamfile:
    full_path = line.strip()
    path, file = os.path.split(full_path)
    if path not in paths:
        paths.append(path)

print len(paths)
print len(files)
print count_files
output_path = 'merged_bams'

for path in paths:
    print 'Agora sim', path
    full_path = '%s' % (path[2:])#.replace(' ', '\ ')

    #os.system('ls -lah %s' % (full_path))
    bam_files = []
    for file in os.listdir(full_path):
        if file.endswith(".bam"):
            #print(file)
            bam_files.append(file)
    print len(bam_files)

    first_bam = bam_files[0]
    print first_bam
    tag = "_".join(first_bam.split('_')[:-1])

    new_bam_files = []
    for count in range(1,len(bam_files)+2):
        #print count
        if count != 23:
            bam_file = full_path.replace(' ', '\ ')+'/'+tag+'_chr%s.hg19.bam' % (count)
            new_bam_files.append(bam_file)

    # print new_bam_files
    command = "samtools merge -f %s.bam %s" % (tag, " ".join(new_bam_files))
    #print command
    os.system(command)

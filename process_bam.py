import os
import hashlib
import pysam
def sizeof_fmt(num, suffix='B'):
    for unit in ['','Ki','Mi','Gi','Ti','Pi','Ei','Zi']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Yi', suffix)


#bams_file = open('/media/raony/storage/lista_bams.txt')
path = '/media/raony/storage/'
#for line in bams_file:
#	print line
#	print sizeof_fmt(os.path.getsize(path+line.strip()))
#	print os.stat(path+line.strip()).st_size 

def md5(fname):
    hash = hashlib.md5()
    with open(fname) as f:
        for chunk in iter(lambda: f.read(4096), ""):
            hash.update(chunk)
    return hash.hexdigest()

bam_list = []
bam_md5 = []
ampliseq_list = []
for root, dirs, files in os.walk(path):
    for name in files:
        # print(os.path.join(root, name))
        if name.endswith('.bam'):
        	bamfile = os.path.join(root, name)
        	bam_list.append(bamfile)
        	#now print file size
        	print bamfile
                print "TAMANHO:"+sizeof_fmt(os.path.getsize(bamfile))
        	#find duplicated
        	#print md5(bamfile)
		samfile = pysam.AlignmentFile(bamfile, "rb")
                print 'ID:'+samfile.header['RG'][0]['ID']              
                if 'IonXpress' in samfile.header['RG'][0]['ID']:
                    print 'achou ampliseq'
                    ampliseq_list.append(bamfile)
    # for name in dirs:
        # print(os.path.join(root, name))
print 'ampliseqs', len(ampliseq_list)

print 'apenas ampliseqs'


def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]

print chunks(ampliseq_list, 3)

chunk_count = 1
for chunk in chunks(ampliseq_list, 15):
    ampliseq_filelist = open('ampliseq.group%s.list' % (chunk_count), 'w')
    chunk_count += 1
    for bam in chunk:
        ampliseq_filelist.writelines(bam+'\n')
    ampliseq_filelist.close()

#for file in ampliseq_list:
#    print file

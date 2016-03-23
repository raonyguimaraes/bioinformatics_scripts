pedigree = open('/projects/1000genomes/integrated_call_sets/integrated_call_samples.20101123.ped', 'r')
header = pedigree.next()
individual_dict = {}
for line in pedigree:
	# print line
	individual = line.split('\t')
	# print individual
	individual_id = individual[1]
	individual_dict[individual_id] = {}
	individual_dict[individual_id]['family'] = individual[0]
	individual_dict[individual_id]['father'] = individual[2]
	individual_dict[individual_id]['mother'] = individual[3]
	individual_dict[individual_id]['country'] = individual[6]

# print individual_dict

# #editped_1092Exomes
pedfile = open('/projects/1000genomes/integrated_call_sets/1092exomes/1092exomes_sureselect.ped', 'r')
outpedfile = open('/projects/1000genomes/integrated_call_sets/1092exomes/1092exomes_sureselect.allinfo.ped', 'w')
for line in pedfile:
 	row = line.split('\t')
 	individual_id = row[1]
 	#edit row data
 	row[0] = individual_dict[individual_id]['family']
 	row[1] = "%s_%s" % (individual_id, individual_dict[individual_id]['country'])
 	if individual_dict[individual_id]['father'] != '0':
 		row[2] = "%s_%s" % (individual_dict[individual_id]['father'], individual_dict[individual_id]['country'])
	if individual_dict[individual_id]['mother'] != '0':
 		row[3] = "%s_%s" % (individual_dict[individual_id]['mother'], individual_dict[individual_id]['country'])
 	#output new file
 	outpedfile.writelines("\t".join(row))
 	
 	
outpedfile.close()

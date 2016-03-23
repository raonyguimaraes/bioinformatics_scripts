#Script to put the names on REAP Result
# Family_id
# Country
# Continent

AMR = ['MXL', 'CLM', 'PUR'] 
AFR = ['ASW', 'LWK', 'YRI']
EAS = ['JPT', 'CHB', 'CHS']
EUR = ['TSI', 'CEU', 'IBS', 'FIN', 'GBR']

continents = {'MXL' : 'AMR', 
 	      'CLM' : 'AMR', 
 	      'PUR' : 'AMR',
 	      'ASW' : 'AFR', 
 	      'LWK' : 'AFR',
 	      'YRI' : 'AFR',
 	      'JPT' : 'EAS',
 	      'CHB' : 'EAS',
 	      'CHS' : 'EAS',
 	      'TSI' : 'EUR',
 	      'CEU' : 'EUR',
 	      'IBS' : 'EUR',
 	      'FIN' : 'EUR',
 	      'GBR' : 'EUR'
 	      }

k1000pedfile = ''
reap_result_file = '/home/raony/mendel/projects/relatedness/bra_trios/triobra5/REAP_pairs_relatedness.txt'


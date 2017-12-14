import os

delfiles = ['alphaPsiVariations', 'alphaPsiVariationsTable', 'scenarios', 'scenariosDescriptive']

for i in range(11):
	for j in range(11):
		delfiles.append(f'psi{i}_alpha{j}_profits')
	
for delfile in delfiles:
	for ext in ['.csv', '.json']:
		print(f'Removing: {delfile+ext}')
		os.remove(delfile+ext)
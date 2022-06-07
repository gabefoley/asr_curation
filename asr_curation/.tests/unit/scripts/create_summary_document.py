import os
import glob
import shutil
base_path = "/".join(snakemake.input.pdf.split("/")[0:-1])

subsets = [x for x in glob.glob(f'{base_path}/{snakemake.wildcards.dataset}/subsets/**/*cleaned.ipynb')]


with open(snakemake.output.config, "w+") as config:
	config.write(f'title: "{snakemake.wildcards.dataset}"\n')
	# config.write(f'logo: images/logo-wide.svg')




with open(snakemake.output.toc, "w+") as toc:
	toc.write(f'format: jb-book \n')
	toc.write(f'root: {snakemake.wildcards.dataset}_summary_cleaned.ipynb \n')
	if subsets:
		toc.write('chapters:\n')
		for subset in subsets:
			toc.write(f'- file: {subset.split("/")[-1]}\n')
	# toc.write(f'chapters:\n')
	# toc.write(f'- file: {snakemake.wildcards.dataset}_summary')

print ('kangaroo')
curr_path = "_".join(snakemake.output.config.split("_")[0:-1])
print (curr_path)


for subset in subsets:
	shutil.copy(subset, curr_path)

shutil.copy(snakemake.input.pdf, curr_path)


print (f'{base_path}/{snakemake.wildcards.dataset}/subsets/**/*cleaned.ipynb')
print ('subsets are ')
print (subsets)
print (os.path.basename(x).split('.')[0] for x in glob.glob(snakemake.output.config.split("_")[0]) )
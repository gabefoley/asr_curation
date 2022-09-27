#!/usr/bin/env python
# coding: utf-8

# In[1]:


######## snakemake preamble start (automatically inserted, do not edit) ########
import sys; sys.path.extend(['/Users/uqgfoley/opt/miniconda3/envs/asr_curation/lib/python3.9/site-packages', '/Users/uqgfoley/Dropbox/Code/Python_Workspace/asr_curation/notebooks']); import pickle; snakemake = pickle.loads(b"\x80\x04\x95\xe8\x07\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8csworkflows/example_workflow/datasets/kari_example_ec_1_1_1_86/subsets/all_kari/kari_example_ec_1_1_1_86_all_kari.aln\x94\x8c\x81workflows/example_workflow/datasets/kari_example_ec_1_1_1_86/subsets/all_kari/csv/kari_example_ec_1_1_1_86_all_kari_alignment.csv\x94e}\x94(\x8c\x06_names\x94}\x94(\x8c\x03aln\x94K\x00N\x86\x94\x8c\x03csv\x94K\x01N\x86\x94u\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x15\x8c\tfunctools\x94\x8c\x07partial\x94\x93\x94h\x06\x8c\x19Namedlist._used_attribute\x94\x93\x94\x85\x94R\x94(h\x1b)}\x94\x8c\x05_name\x94h\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh\x0fh\nh\x11h\x0bub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94\x8c\x99workflows/example_workflow/datasets/kari_example_ec_1_1_1_86/dataset_summary/kari_example_ec_1_1_1_86/subsets/all_kari/temp/all_kari_subset_summary.ipynb\x94a}\x94(h\r}\x94\x8c\x07summary\x94K\x00N\x86\x94sh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh,h)ub\x8c\x06params\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94]\x94(\x8c\x14lineage_superkingdom\x94\x8c\x0bxref_supfam\x94\x8c\x0cxref_panther\x94\x8c\nKARI_Class\x94\x8c\x02ec\x94ea}\x94(h\r}\x94\x8c\x0fannotation_cols\x94K\x00N\x86\x94sh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bhCh;ub\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94(\x8c\x18kari_example_ec_1_1_1_86\x94\x8c\x08all_kari\x94e}\x94(h\r}\x94(\x8c\x07dataset\x94K\x00N\x86\x94\x8c\x06subset\x94K\x01N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94b\x8c\x07dataset\x94hR\x8c\x06subset\x94hSub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c0/var/folders/xs/24s9hwqd191f2x7rhdy6_ryc0000gr/T\x94e}\x94(h\r}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bhmK\x01hoK\x01hqhjub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94\x8c\x99workflows/example_workflow/datasets/kari_example_ec_1_1_1_86/dataset_summary/kari_example_ec_1_1_1_86/subsets/all_kari/temp/all_kari_subset_summary.ipynb\x94a}\x94(h\r}\x94\x8c\x08notebook\x94K\x00N\x86\x94sh\x13]\x94(h\x15h\x16eh\x15h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x15sNt\x94bh\x16h\x19h\x1b\x85\x94R\x94(h\x1b)}\x94h\x1fh\x16sNt\x94bh\x83h\x80ub\x8c\x06config\x94}\x94(\x8c\x07workdir\x94\x8c\x1aworkflows/example_workflow\x94\x8c\x08fastadir\x94\x8c workflows/example_workflow/fasta\x94\x8c\x06subdir\x94\x8c'workflows/example_workflow/subset_rules\x94\x8c\x0fannotation_cols\x94]\x94(h<h=h>h?h@e\x8c\x10blocked_datasets\x94]\x94u\x8c\x04rule\x94\x8c\x15create_subset_summary\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cD/Users/uqgfoley/Dropbox/Code/Python_Workspace/asr_curation/notebooks\x94ub."); from snakemake.logging import logger; logger.printshellcmds = False; import os; os.chdir(r'/Users/uqgfoley/Dropbox/Code/Python_Workspace/asr_curation');
######## snakemake preamble end #########


# In[2]:


import matplotlib.pyplot as plt
from collections import defaultdict
import py3Dmol
import pandas as pd
import re
import requests
from IPython.display import display, Markdown


# In[3]:


base_url = "https://www.ebi.ac.uk/pdbe/"

api_base = base_url + "api/"

summary_url = api_base + 'pdb/entry/summary/'
secondary_structure_url = api_base + 'pdb/entry/secondary_structure/'
ligand_url = api_base + '/pdb/entry/ligand_monomers/'


def make_summary(data):
    """
    This function creates a summary for a PDB entry
    by getting data for an entry, and extracting
    pieces of information
    
    :param data: Dict
    :return: String
    """
    
    pdb_id = ""
    # Certain calls could return multiple PDB entries,
    # but the GET summary call we use in this exercise
    # will always return only one PDB entry
    for key in data.keys():
        pdb_id = key

    # The data is a list of dictionaries, and for the summary information,
    # it is always the first element of the list
    entry = data[pdb_id][0]
    
    # Getting the title of the entry
    title = entry['title']
    
    # Getting the release date of the entry
    release_date = entry['release_date']
    # Formatting the entry to make it more user-friendly
    formatted_release_date = "%s/%s/%s" % (
        release_date[:4], 
        release_date[4:6], 
        release_date[6:])
    
    # Getting the experimental methods
    # Note that there can be multiple methods, so this is a list that
    # needs to be iterated
    experimental_methods = ""
    for experimental_method in entry["experimental_method"]:
        if experimental_methods:
            experimental_methods += " and "
        experimental_methods += experimental_method
        
        
        
    # Getting the assemblies
    assemblies = ''
    for assembly in entry['assemblies']:
        if assembly:
            # Blank out the form if it is a monomer (because it must be homo form)
            if assembly['name'] == 'monomer':
                form = ""
            else:
                form = assembly['form']
                
            if assembly['preferred']:
                
                assemblies += f"Preferred form with assembly ID {assembly['assembly_id']} is a {form}{assembly['name']}\n"
            else:
                assemblies += f"Non-preferred form with assembly ID {assembly['assembly_id']} is a {assembly['form']}{assembly['name']}\n"

    # Getting the author list
    authors = "".join(entry['entry_authors'])
    
#         'assemblies': [{'assembly_id': '1', 'form': 'homo', 'preferred': True, 'name': 'monomer'}, {'assembly_id': '2', 'form': 'homo', 'preferred': False, 'name': 'dimer'}]}]}
    
    
        
     
    # Creating the summary text using all the extracted 
    # information
    summary = f'Entry is titled "{title}" and was released on {formatted_release_date} \n\nThis entry was determined using {experimental_methods} \nThe authors are {authors} \n\nAssembly information -\n{assemblies}'

    return summary

def make_request(url, mode, pdb_id):
    """
    This function can make GET and POST requests to
    the PDBe API
    
    :param url: String,
    :param mode: String,
    :param pdb_id: String
    :return: JSON or None
    """
    if mode == "get":
        response = requests.get(url=url+pdb_id)
    elif mode == "post":
        response = requests.post(url, data=pdb_id)

    if response.status_code == 200:
        return response.json()
    else:
        print("[No data retrieved - %s] %s" % (response.status_code, response.text))
    
    return None

def get_secondary_structure_ranges(pdb_id=None, pdb_list=None):
    """
    This function calls the PDBe API and retrieves the residue
    ranges of secondary structural elements in a single PDB entry
    or in a list of PDB entries
    
    :param pdb_id: String,
    :param pdb_list: String
    :return: None
    """
    # If neither a single PDB id, nor a list was provided,
    # exit the function
    if not pdb_id and not pdb_list:
        print("Either provide one PDB id, or a list of ids")
        return None
    
    if pdb_id:
        # If a single PDB id was provided, call the API with GET
        data = make_request(secondary_structure_url, "get", pdb_id)
    else:
        # If multiple PDB ids were provided, call the API with POST
        # The POST API call expects PDB ids as a comma-separated lise
        pdb_list_string = ", ".join(pdb_list)
        data = make_request(secondary_structure_url, "post", pdb_list_string)
        
    # When no data is returned by the API, exit the function
    if not data:
        print("No data available")
        return None
    
    # Loop through all the PDB entries in the retrieved data
    for entry_id in data.keys():
        entry = data[entry_id]
        molecules = entry["molecules"]
        
        # Loop through all the molecules of a given PDB entry
        for i in range(len(molecules)):
            chains = molecules[i]["chains"]
            
            # Loop through all the chains of a given molecules
            for j in range(len(chains)):
                secondary_structure = chains[j]["secondary_structure"]
                helices = secondary_structure["helices"] if 'helices' in secondary_structure else []
                strands = secondary_structure["strands"] if 'strands' in secondary_structure else []
                helix_list = []
                strand_list = []
                
                # Loop through all the helices of a given chain
                for k in range(len(helices)):
                    start = helices[k]["start"]["residue_number"]
                    end = helices[k]["end"]["residue_number"]
                    helix_list.append("%s-%s" % (start, end))
                
                # Loop through all the strands of a given chain
                for l in range(len(strands)):
                    start = strands[l]["start"]["residue_number"]
                    end = strands[l]["end"]["residue_number"]
                    strand_list.append("%s-%s" % (start, end))
                    
                report = "%s chain %s has " % (entry_id, chains[j]["chain_id"])
                if len(helix_list) > 0:
                    report += "helices at residue ranges %s " % str(helix_list)
                else:
                    report += "no helices "
                report += "and "
                if len(strand_list) > 0:
                    report += "strands at %s" % str(strand_list)
                else:
                    "no strands"
                print(report)
                
    return None

def get_ligand_information(data):
    pdb_id = ""
    
    summary = ""

    for key in data.keys():
        pdb_id = key
        
    for ligand in data[pdb_id]:
        summary += f"Ligand : {ligand['chem_comp_name']} at position {ligand['author_residue_number']}\n"
        
    return summary
    

def get_entry_from_api(pdb_id, api_url):
    """
    This function will make a call to the PDBe API using
    the PDB id and API url provided as arguments
    
    :param pdb_id: String,
    :param api_url: String
    :return: Dict or None
    """
    if not re.match("[0-9][A-Za-z][A-Za-z0-9]{2}", pdb_id):
        print("Invalid PDB id")
        return None
    
    # Make a GET call to the API URL
    get_request = requests.get(url=api_url+pdb_id)
    
    if get_request.status_code == 200:
        # If there is data returned (with HTML status code 200)
        # then return the data in JSON format
        return get_request.json()
    else:
        # If there is no data, print status code and response
        print(get_request.status_code, get_request.text)
        return None


# As you can hopefully see, the data displayed is very similar to
# what we had in the mock data in previous sections - however,
# this is actual data coming from the PDBe API
    


# In[4]:


print(display(Markdown(f'# Subset - {snakemake.wildcards.subset}')))


# In[5]:


print (f'Subset sheet - {snakemake.wildcards.subset}')


# In[6]:


df = pd.read_csv(snakemake.input.csv)
non_fragment_df = df[df['fragment'] == False]
entry_df = df.dropna(axis=1, how='all')


# # EC numbers

# In[7]:


ec_nums = list(pd.unique(df['ec']))

ec_set = set()

for num in ec_nums:
    for split in str(num).split(";"):
        ec_set.add(split.strip())
    
print (f'The EC numbers found in the current data set are {ec_set}')


# # Taxonomic distribution

# In[8]:


taxonomy_cols = ['lineage_superkingdom', 'lineage_phylum']

for col in taxonomy_cols:
    fig, ax = plt.subplots(figsize=(len(col) / 4 ,10))
    chart = df[col].value_counts().plot.barh(title=col, ax=ax)
    plt.show()


# # Annotation distributions

# In[9]:


key_annotation_cols = snakemake.params.annotation_cols

for col in key_annotation_cols:
    if col in entry_df:
        fig, ax = plt.subplots(figsize=(len(col) / 3 ,10))
        chart = entry_df[col].value_counts().plot.barh(title=col, ax=ax)
        plt.show()


# In[10]:


seq_len = df['length'].describe()
non_fragment_seq_len = non_fragment_df['length'].describe()



def print_sequence_summarys(seq_summary):
    print (f"Number of sequences : {seq_summary['count']}")
    print (f"Smallest sequnce length : {seq_summary['min']}")
    print (f"Longest sequnce length : {seq_summary['max']}")
    print (f"Average sequnce length : {seq_summary['mean']}")

print ("Sequence statistics for all sequences")
print_sequence_summarys(seq_len)

print ()

if non_fragment_seq_len.all():

    print ("Sequence statistics for sequences (non-fragments)")
    print_sequence_summarys(non_fragment_seq_len)


# # Experimental data (BRENDA)

# In[11]:


skip_brenda_cols = [
'BRENDA_CL', 'BRENDA_GI', 
'BRENDA_SN', 'BRENDA_SY', 
'BRENDA_MW', 'BRENDA_SP',
'BRENDA_NSP', 'BRENDA_PM',
'BRENDA_LO', 'BRENDA_SU', 
'BRENDA_PU', 'BRENDA_ST',
'BRENDA_CR', 'BRENDA_CF',
'BRENDA_RN', 'BRENDA_RT',
'BRENDA_ME', 'BRENDA_REFERENCES']


# Only get the columns that start with BRENDA and that we don't want to skip
brenda_cols =[x for x in df if x.startswith("BRENDA") and 'COMMENT' not in x and 'REFS' not in x and 'UNITS' not in x and not any(skip in x for skip in skip_brenda_cols)]

# Drop rows if they don't have at least one entry in one of the BRENDA columns
b_df = df[brenda_cols].dropna(thresh=1)
print (len(brenda_cols))
fig, ax = plt.subplots(figsize=(len(brenda_cols) / 4 ,10))


# Create and save a plot of the BRENDA column counts
if not b_df.empty:
    b_df.count().plot.barh(ax=ax)
    display()


# # Structural information

# In[12]:


pdb_df = df[['accession', 'xref_pdb']].dropna()
pdb_dict = defaultdict(list)
for pdb_entries in pdb_df.itertuples():
    pdb_dict[pdb_entries.accession] = [x for x in pdb_entries.xref_pdb.strip().split(";") if len(x) > 0]

total_structures = sum([len(x) for x in pdb_dict.values()])

print (f"There are {total_structures} total structures across {len(pdb_dict.keys())} unique sequences")

print ("Looking at the first few...")

sample_pdb_dict = {k: pdb_dict[k] for k in list(pdb_dict)[:3]}

for prot_id, pdb_ids in sample_pdb_dict.items():
    for pdb_id in pdb_ids:
        
        print(display(Markdown(f'# {pdb_id} PDB entry information')))

        print(make_summary(get_entry_from_api(pdb_id, summary_url)))
        
        display(Markdown(f'# {pdb_id} Ligand informtion'))    
        print(get_ligand_information(get_entry_from_api(pdb_id, ligand_url)))
        
    
        display(Markdown(f'# {pdb_id} Secondary structure ranges'))
        get_secondary_structure_ranges(pdb_id)
        
        display(Markdown(f'# {pdb_id}'))
        
        view = py3Dmol.view(query=f'pdb:{pdb_id}')
        view.setStyle({'cartoon':{'color':'spectrum'}})
        display(view)


# <!-- # Alignment information -->

# In[13]:


# Alignment information (not working at the moment)
# print ('Make this work on server')
# import sequence
# subset_aln = sequence.readFastaFile(f"{snakemake.input.aln}")
# display(Markdown(f'Alignment has a length of {len(subset_aln[0].sequence)} positions'))


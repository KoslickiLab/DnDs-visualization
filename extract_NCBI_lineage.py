from Bio import Entrez
from tqdm import tqdm
import time
import os
import pandas as pd



# Set your email (NCBI requires this)
Entrez.email = "your.email@example.com"

# Define the desired ranks
desired_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def get_lineage_and_ranks(assembly_id):
    try:
        # Search for the assembly in NCBI
        search_handle = Entrez.esearch(db="assembly", term=assembly_id)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        if not search_results['IdList']:
            return None, None

        # Get the assembly UID
        assembly_uid = search_results['IdList'][0]

        # Fetch the assembly summary to get the taxonomy ID
        summary_handle = Entrez.esummary(db="assembly", id=assembly_uid)
        summary = Entrez.read(summary_handle)
        summary_handle.close()

        tax_id = summary['DocumentSummarySet']['DocumentSummary'][0]['Taxid']

        # Fetch the taxonomy information
        tax_handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
        tax_info = Entrez.read(tax_handle)
        tax_handle.close()

        # Extract the lineage and ranks
        lineage = [lineage['ScientificName'] for lineage in tax_info[0]['LineageEx']]
        ranks = [lineage['Rank'] for lineage in tax_info[0]['LineageEx']]
        
        # Add the current taxon
        lineage.append(tax_info[0]['ScientificName'])
        ranks.append(tax_info[0]['Rank'])
        
        # Filter the lineage and ranks to only include desired ranks
        filtered_lineage = []
        filtered_ranks = []
        for rank, name in zip(ranks, lineage):
            if rank in desired_ranks:
                filtered_lineage.append(name)
                filtered_ranks.append(rank)

        # Introduce a delay to avoid hitting rate limits
        time.sleep(1)

        return filtered_lineage, filtered_ranks

    except Exception as e:
        if "429" in str(e):
            print(f"Rate limit exceeded, retrying after a delay... ({assembly_id})")
            time.sleep(10)  # Wait 10 seconds before retrying
            return get_lineage_and_ranks(assembly_id)
        else:
            print(f"An error occurred: {e}")
            return None, None

# Read Comparison File
file_path = '/home/other/jpf5265/DnDs-visualization/full_dataset/fmh_omega_7.csv'
file = pd.read_csv(file_path, sep=',', header=0)
query_name_list = ['_'.join(x.split('_')[1:]) for x in file['query_name_x'].to_list()]
match_names_list = ['_'.join(x.split('_')[1:]) for x in file['match_name_x'].to_list()] 
assembly_ids = list(set(query_name_list + match_names_list))

# set up a flag to check if there is a difference between two runs
has_diff = 1

# Read lineage info result if it exists
if os.path.exists("lineage_information.csv"):
    old_df = pd.read_csv("lineage_information.csv")
    existing_assembly_ids = set(old_df['Assembly ID'].to_list())
    run_assembly_ids = list(set(assembly_ids).difference(existing_assembly_ids))
else:
    run_assembly_ids = assembly_ids
    
# set up some checking variables
last_diff_num = len(run_assembly_ids)
iter = 1  

while has_diff == 1 and last_diff_num != 0:

    print(f"Running iteration {iter}...", flush=True)
    
    # Get lineage info
    lineages_ranks = {assembly_id: get_lineage_and_ranks(assembly_id) for assembly_id in tqdm(run_assembly_ids)}

    # Filter out entries with no results
    lineages_ranks = {k: v for k, v in lineages_ranks.items() if v[0] is not None}

    # Create an empty DataFrame with N rows and 8 columns
    new_df = pd.DataFrame(index=range(len(lineages_ranks)), columns=range(8))

    # add column names
    new_df.columns = ['Assembly ID'] + desired_ranks

    # add lineage information into the dataframe
    for i, (assembly_id, (lineage, ranks)) in enumerate(lineages_ranks.items()):
        new_df.loc[i, 'Assembly ID'] = assembly_id
        for j, (rank, name) in enumerate(zip(ranks, lineage)):
            if rank in desired_ranks:
                new_df.loc[i, rank] = name

    # Add additional column "root" which is always "cellular organisms", betwen “Assembly ID” and “superkingdom”
    new_df.insert(1, "root", "cellular organisms")

    # Display the DataFrame
    print(new_df)

    # Combine the old and new dataframe
    if 'old_df' not in locals():
        df = new_df
    else:
        df = pd.concat([old_df, new_df], ignore_index=True)

    # Check difference again
    existing_assembly_ids = set(df['Assembly ID'].to_list())
    run_assembly_ids = list(set(assembly_ids).difference(existing_assembly_ids))
    this_diff_num = len(run_assembly_ids)
    if this_diff_num == last_diff_num:
        has_diff = 0
        # Save to CSV
        df.to_csv("lineage_information.csv", index=False)
        iter = 1
        break
    else:
        last_diff_num = this_diff_num
        old_df = df
        # add iter by 1
        iter += 1

if iter != 1:
    # save to CSV
    df.to_csv("lineage_information.csv", index=False)   

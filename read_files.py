import pandas as pd
import numpy as np
import os


def read_file(protein_file="", peptide_file="", file_id=0):
    """ Process and populate pandas DataFrames with the data from 
    the protein and peptide files, saving only the Protein.Names, 
    Precursor.Id (for peptide), Protein.Group (for protein), and 
    the file name columns. Also rename the run names with a run ID,
    and save this infor

    Parameters:
        protein_file (string): A pandas readable file with the protein data
        peptide_file (string): A pandas readable file with the peptide data
        file_id (int): The number denoting the group in the run ID

    Returns:
        data_obj (dict): A dictionary with the following keys/values
            "run_metadata": A pandas DataFrame mapping run names to ids
            "pep_abundance": A pandas DataFrame with peptide data
            "prot_abundance": A pandas DataFrame with protein data
    """
    ## Create a pandas DataFrame for peptide and protein ##
    try:
        protein_table = pd.read_table(protein_file, low_memory=False)
    except FileNotFoundError:
        return {"early": f"Protein file{file_id} '{protein_file}' does not exist."}

    try:
        peptide_table = pd.read_table(peptide_file, low_memory=False)
    except FileNotFoundError:
        return {"early": f"Peptide file{file_id} '{peptide_file}' does not exist."}
    
    # Filters out proteins with "contam" or that are null
    peptide_table = peptide_table[~peptide_table['Protein.Group'].str.contains("contam", na=False)]
    protein_table = protein_table[~protein_table['Protein.Group'].str.contains("contam", na=False)]

    # Separate the columns of just the file names, 
    # "Precursor.Id" for peptides, "Protein.Group" for protein (these are unique), 
    # and "Protein.Names" (used to filter by organism)
    peptide_cols = peptide_table.filter(regex='\\\\|Precursor.Id|Protein.Names').columns
    protein_cols = protein_table.filter(regex="\\\\|Protein.Group|Protein.Names").columns

    # Create the abundance matrix
    pep_abundance = peptide_table.loc[:, peptide_cols]
    prot_abundance = protein_table.loc[:, protein_cols]


    ## Create the run IDs ##
    # Create a list of just the full file path names
    full_pep_run_names = peptide_cols.drop(['Precursor.Id', 'Protein.Names']).to_list()
    full_prot_run_names = protein_cols.drop(['Protein.Group', 'Protein.Names']).to_list()

    # If .raw is in the file name, strip off everything after .raw.
    if ".raw" in full_pep_run_names[0]:
        run_names = [name.split(".raw")[0] for name in full_pep_run_names]
    else:
        run_names = [os.path.splitext(name)[0] for name in full_pep_run_names]

    # Create a list of the run_IDs
    num_files = len(run_names)
    ID_generator = lambda x: str(file_id) + "-" + str(x)
    run_IDs = [ID_generator(i) for i in range(num_files)]

    # Create the table mapping run_name to run_ID
    run_metadata = pd.DataFrame({"run_name": np.array(run_names), 
                                "run_ID": np.array(run_IDs)})
    

    ## Rename the DataFrames with run_IDs ##
    # Create a dictionary key run_name and value run_ID
    pep_rename_dict = dict(zip(full_pep_run_names, run_IDs))
    prot_rename_dict = dict(zip(full_prot_run_names, run_IDs))

    pep_abundance.rename(columns=pep_rename_dict, inplace=True)
    prot_abundance.rename(columns=prot_rename_dict, inplace=True)
    
    return {"run_metadata": run_metadata,
            "pep_abundance": pep_abundance, 
            "prot_abundance": prot_abundance}



def read_files(filelist = []):
    """ Read in the combined data from multiple runs. Call read_file 
    to process and populate pandas DataFrames with the data, saving 
    only the Protein.Names, Precursor.Id (for peptide), Protein.Group 
    (for protein), and the file name columns. Also call group_data to 
    combine the data into one data object.

    Parameters:
        filelist (list): A list of dictionaries that contain the protein
            and peptide matrices to analyze

    Returns:
        data_obj (dict): A dictionary with the following keys/values
            "run_metadata": A pandas DataFrame mapping run names to ids
            "pep_abundance": A pandas DataFrame with peptide data
            "prot_abundance": A pandas DataFrame with protein data
    """
    # Initiate the data_objects list
    data_objects = []

    # Initiate i as a file_id counter
    i = 0

    for file_dict in filelist:
        # Give error messages if syntax is not correct
        if "protein_file" not in file_dict:
            print('Dictionary Syntax incorrect. Must include "protein_file".')
            return
        if "peptide_file" not in file_dict:
            print('Dictionary Syntax incorrect. Must include "peptide_file".')
            return
        
        protein_file = file_dict["protein_file"]
        peptide_file = file_dict["peptide_file"]

        # Read the file and append it to data_objects
        data_obj = read_file(protein_file, peptide_file, i)

        # Give error message if returned early
        if "early" in data_obj:
            print(data_obj["early"])
            return

        data_objects.append(data_obj)
        i += 1
        



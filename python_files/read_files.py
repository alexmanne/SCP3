import pandas as pd
import numpy as np
import os
import yaml
from time import time



def diann_1_9_tables(file="", is_protein=True):
    """ Process and populate pandas DataFrames with the data from 
    the DIA-NN 1.9 files (protein or peptide), saving only certain columns.

    **Protein Columns**: Protein.Names, Protein.Group, file names     
    **Peptide Columns**: Protein.Names, Precursor.Id, file names     
   
    Also rename these columns.     
    - Protein.Names: *Protein Name*
    - Protein.Group: *Accession*
    - Precursor.Id: *Sequence*

    Parameters:
        file (string):  A pandas readable file with the data      
        is_protein (bool):  True if protein file, False if peptide

    Returns:
        abundance (DataFrame):  pandas DataFrame with filtered data 

    Raises:
        FileNotFoundError:  If file does not exist
    """
    ## Verify that the files are valid ##
    try:
        df_table = pd.read_table(file, low_memory=False)
    except FileNotFoundError:
        print(f"File '{file}' does not exist.")
        raise FileNotFoundError(f"File '{file}' does not exist.")
    

    ## Create the base abundance table ##
    # Filters out proteins with "contam" or that are null
    df_table = df_table[~df_table['Protein.Group'].str.contains("contam", na=False)]

    # Separate the columns of just the file names, 
    # "Precursor.Id" for peptides, "Protein.Group" for protein (these are unique), 
    # and "Protein.Names" (used to filter by organism)
    if is_protein:
        col_names = df_table.filter(regex="\\\\|Protein.Group|Protein.Names").columns
        identifier = "Protein.Group"
        new_name = "Accession"
    else:
        col_names = df_table.filter(regex='\\\\|Precursor.Id|Protein.Names').columns
        identifier = "Precursor.Id"
        new_name = "Sequence"

    # Create the abundance matrix
    abundance = df_table.loc[:, col_names]


    ## Rename the columns ##
    # Create a list of just the full file path names
    run_names = col_names.drop([identifier, 'Protein.Names']).to_list()

    # If .raw is in the file name, strip off everything after .raw.
    if ".raw" in run_names[0]:
        short_run_names = [name.split(".raw")[0] for name in run_names]
    else:
        short_run_names = [os.path.splitext(name)[0] for name in run_names]

    # Create the renaming dictionary
    rename_dict = dict(zip(run_names, short_run_names))
    rename_dict[identifier] = new_name
    rename_dict["Protein.Names"] = "Protein Name"

    # Rename the file columns
    abundance.rename(columns=rename_dict, inplace=True)

    return abundance


def fragpipe_22_tables(file="", is_protein=True, 
                       min_unique_peptides=1, use_maxlfq=False):
    """ Process and populate pandas DataFrames with the data from the 
    FragPipe v22.0 files (protein and ion), saving only certain columns.
    
    **Protein Columns**: Entry Name, Protein ID, file names     
    **Peptide Columns**: Entry Name, Modified Sequence, file names   

    Also rename these columns.     
    - Entry Name: *Protein Name*
    - Protein ID: *Accession*
    - Modified Sequence: *Sequence*

    Parameters:
        file (string): A pandas readable file with the data
        is_protein (bool): True if protein file, False if peptide
        min_unique_peptides (int): The minimum number of unique peptides
        use_maxlfq (bool): Determines whether or not to use MaxLFQ Intensity     
            Default: False (uses "Intensity" columns without MaxLFQ)     
            Only supported on Protein files

    Returns:
        abundance (DataFrame): pandas DataFrame with the filtered data

    Raises:
        FileNotFoundError: If file does not exist     
                           If not the combined_ion.tsv file
    """
    ## Verify that the files are valid ##
    try:
        df_table = pd.read_table(file, low_memory=False)
    except FileNotFoundError:
        print(f"File '{file}' does not exist.")
        raise FileNotFoundError(f"File '{file}' does not exist.")

    if not is_protein and "M/Z" not in df_table:
        print(f"Peptide file must be the combined_ion.tsv file.")
        raise FileNotFoundError(f"Peptide file must be the combined_ion.tsv file.")
    

    ## Create the base abundance table ##
    # Filters out proteins with "contam" or that are null
    df_table = df_table[~df_table['Protein'].str.contains("contam", na=False)]

    # Separate the columns of
    # "Peptide Sequence" (unique) for peptides, 
    # "Protein ID" (unique) for protein, and
    # "Entry Name" (used to filter by organism)

    # Edit the table and save the protein specific values
    if is_protein:

        # Filter out the proteins with total peptide count less than min_unique_peptides
        df_table.query("`Combined Total Peptides` >= @min_unique_peptides", inplace=True)

        # Get the column names
        col_names = df_table.filter(regex="Protein ID|Entry Name").columns.to_list()

        # Add the Intensity columns (depending on use_maxlfq)
        if use_maxlfq is True:
            col_names += (df_table.filter(regex=' MaxLFQ Intensity').columns.to_list())
        else:
            col_names += (df_table.filter(regex='(?<!MaxLFQ) Intensity').columns.to_list())

        # Save these values for later
        identifier = "Protein ID"
        new_name = "Accession"

    # Edit the table and save the peptide specific values
    else:
        # Add the Charge number to the peptide sequence
        df_table["Modified Sequence"] = df_table["Modified Sequence"] + df_table["Charge"].astype(str)
        col_names = df_table.filter(regex='Modified Sequence|Entry Name|(?<!MaxLFQ) Intensity').columns.to_list()

        # Save these values for later
        identifier = 'Modified Sequence'
        new_name = "Sequence"

    # Create the abundance data frame
    abundance = df_table.loc[:, col_names]


    # ## Rename the columns ##
    # Remove the column names that are not run names 
    col_names.remove(identifier)
    col_names.remove('Entry Name')

    # Get just the run name (Assumes the run name has no spaces)
    run_names = [name.split(" ")[0] for name in col_names]
    
    # Create the renaming dictionary
    rename_dict = dict(zip(col_names, run_names))
    rename_dict[identifier] = new_name
    rename_dict["Entry Name"] = "Protein Name"

    # Rename the file columns
    abundance.rename(columns=rename_dict, inplace=True)

    return abundance


def read_file(protein_file="", peptide_file="", file_id=0, 
              processing_app='', min_unique_peptides=1, use_maxlfq=False):
    """ Process and populate pandas DataFrames with the data from 
    the protein and peptide files, saving only the Organism, 
    Sequence (for peptide), Accession (for protein), and the 
    file name columns. Also rename the run names with a run ID,
    and save this information

    Parameters:
        protein_file (string): A pandas readable file with the protein data
        peptide_file (string): A pandas readable file with the peptide data
        file_id (int): The number denoting the group in the run ID
        processing_app (string): Denotes the application used for processing
            Currently supports: "diann", "fragpipe"
        min_unique_peptides (int): The minimum number of unique peptides (used in FragPipe)

    Returns:
        data_obj (dict): A dictionary with the following keys/values
            "run_metadata": A pandas DataFrame mapping run names to ids
            "pep_abundance": A pandas DataFrame with peptide data
            "prot_abundance": A pandas DataFrame with protein data
    
    Raises:
        ValueError: If the processing app is not accepted
        FileNotFoundError: If protein_file or peptide_file does not exist
    """
    ## Verify the processing application ##
    supported_apps = ["diann", "fragpipe"]
    if processing_app not in supported_apps:
        print(f"Processing app {processing_app} not currently supported.")
        raise ValueError(f"Processing app {processing_app} not currently supported.")
    

    ## Create a pandas DataFrame for peptide and protein abundance ##
    prot_abundance = pd.DataFrame()
    pep_abundance = pd.DataFrame()

    if processing_app == "diann":
        if protein_file:
            prot_abundance = diann_1_9_tables(protein_file, is_protein=True)
        if peptide_file:
            pep_abundance = diann_1_9_tables(peptide_file, is_protein=False)
    
    elif processing_app == "fragpipe":
        if protein_file:
            prot_abundance = fragpipe_22_tables(protein_file, is_protein=True, 
                                                min_unique_peptides=min_unique_peptides, 
                                                use_maxlfq=use_maxlfq)
        if peptide_file:
            pep_abundance = fragpipe_22_tables(peptide_file, is_protein=False)

    
    ## Rename the DataFrames with run_IDs ##
    # Create a list of the run_names using pep_abundance by default
    if not pep_abundance.empty:
        run_names = pep_abundance.columns.drop(["Sequence", "Protein Name"]).to_list()
    else:
        run_names = prot_abundance.columns.drop(["Accession", "Protein Name"]).to_list()
    
    num_files = len(run_names)
    ID_generator = lambda x: str(file_id) + "-" + str(x)
    run_IDs = [ID_generator(i) for i in range(num_files)]

    # Create the table mapping run_name to run_ID
    run_metadata = pd.DataFrame({"run_name": np.array(run_names), 
                                "run_ID": np.array(run_IDs)})

    # Create a dictionary key run_name and value run_ID
    rename_dict = dict(zip(run_names, run_IDs))

    pep_abundance.rename(columns=rename_dict, inplace=True)
    prot_abundance.rename(columns=rename_dict, inplace=True)
    
    return {"run_metadata": run_metadata,
            "pep_abundance": pep_abundance, 
            "prot_abundance": prot_abundance}




def group_data(data_objects):
    """ Take the list of data objects and combine them into one.

    Parameters:
        data_objects (list): A list of dictionaries that contain the data 
            objects to combine. Each data object has the following keys/values
                "run_metadata": A pandas DataFrame mapping run names to ids
                "pep_abundance": A pandas DataFrame with peptide data
                "prot_abundance": A pandas DataFrame with protein data

    Returns:
        data_obj (dict): A dictionary with the following keys/values
            "run_metadata": A pandas DataFrame mapping run names to ids
            "pep_abundance": A pandas DataFrame with peptide data
            "prot_abundance": A pandas DataFrame with protein data
    
    Raises:
        KeyError: If one of the keys is missing from a data object
    """

    # Initialize a counter for error catching
    i = 0

    for data_obj in data_objects:

        # Verify the keys are what we want them to be
        for key in data_obj.keys():
            if key not in ["run_metadata", "pep_abundance", "prot_abundance"]:
                print(f"{key} not found in the {i} data object")
                raise KeyError(f"{key} not found in the {i} data object")
            
        # Initialize the final data object if it is the first one. (if i==0)
        if i == 0:
            final_data_obj = data_obj

        # append the information to the final_data_object for the rest
        else:
            # Concatenate the run metadata and reset the index
            final_data_obj["run_metadata"] = pd.concat([final_data_obj["run_metadata"], 
                                                        data_obj["run_metadata"]])
            final_data_obj["run_metadata"].reset_index(drop=True, inplace=True)

            # Outer join the peptide matrices on the Sequence column
            final_data_obj["pep_abundance"] = pd.merge(final_data_obj["pep_abundance"], 
                                                       data_obj["pep_abundance"],
                                                       how="outer", on=["Sequence", "Protein Name"])
            
            # Outer join the protein matrices on the Accession column
            final_data_obj["prot_abundance"] = pd.merge(final_data_obj["prot_abundance"], 
                                                       data_obj["prot_abundance"],
                                                       how="outer", on=["Accession", "Protein Name"])
            
        i += 1

     # Replace 0 values with nan to decrease complexity
    final_data_obj["pep_abundance"].replace(0, np.nan, inplace=True)
    final_data_obj["prot_abundance"].replace(0, np.nan, inplace=True)

    return final_data_obj





def read_files(yaml_file):
    """ Read in the combined data from multiple runs. Call read_file 
    to process and populate pandas DataFrames with the data, saving 
    only the Organism, Peptide Sequence (for peptide), Protein ID 
    (for protein), and the file name columns. Also call group_data to 
    combine the data into one data object.

    Parameters:
        yaml_file (string): A string with the name of the yaml file
            containing the names of peptide/protein files

    Returns:
        data_obj (dict): A dictionary with the following keys/values
            "run_metadata": A pandas DataFrame mapping run names to ids
            "pep_abundance": A pandas DataFrame with peptide data
            "prot_abundance": A pandas DataFrame with protein data
    
    Raises:
        ValueError: If the processing app is not accepted
                    If the dictionary is not correct
        FileNotFoundError: If the yaml_file, protein_file, or 
                    peptide_file do not exist
    """
    ## Verify that the file is valid ##
    try:
        with open(yaml_file, "r") as f:
            filelist = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"Protein file '{yaml_file}' does not exist.")
        raise FileNotFoundError(f"Protein file '{yaml_file}' does not exist.")

    # Initiate the data_objects list
    data_objects = []

    # Initiate i as a file_id counter
    i = 0

    for file_dict in filelist:
        # Give error messages if syntax is not correct
        if "protein_file" not in file_dict:
            print('yaml syntax incorrect. Must include "protein_file".')
            raise ValueError('yaml syntax incorrect. Must include "protein_file".')
        if "peptide_file" not in file_dict:
            print('yaml syntax incorrect. Must include "peptide_file".')
            raise ValueError('yaml syntax incorrect. Must include "peptide_file".')
        if "processing_app" not in file_dict:
            print('yaml syntax incorrect. Must include "processing_app".')
            raise ValueError('yaml syntax incorrect. Must include "processing_app".')
        
    
        ## currently set min peptides to 1. Can be changed with settings
        min_unique_peptides = 5
        ## currently set use_maxlfq to False. Can be changed with settings
        use_maxlfq = False

        # Read the file and append it to data_objects
        data_obj = read_file(protein_file=file_dict["protein_file"], 
                             peptide_file=file_dict["peptide_file"], 
                             file_id=i, 
                             processing_app= file_dict["processing_app"],
                             min_unique_peptides=min_unique_peptides,
                             use_maxlfq=use_maxlfq)
        
        data_objects.append(data_obj)
        i += 1

    # Call group data to put all of the data into one data object
    return group_data(data_objects)
        



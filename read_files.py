import pandas as pd
import numpy as np
import os


def diann_1_9_tables(protein_file="", peptide_file=""):
    """ Process and populate pandas DataFrames with the data from 
    the DIA-NN 1.9 protein and peptide files, saving only the Protein.Names, 
    Precursor.Id (for peptide), Protein.Group (for protein), and the 
    file name columns. Also rename these columns.

    Parameters:
        protein_file (string): A pandas readable file with the protein data
        peptide_file (string): A pandas readable file with the peptide data

    Returns:
        prot_abundance (DataFrame): pandas DataFrame with protein data
            Column names: "Organism", "Accession", and the run names
        pep_abundance (DataFrame): pandas DataFrame with peptide data
            Column names: "Organism", "Sequence", and the run names

    Raises:
        FileNotFoundError: If protein_file or peptide_file does not exist
    """
    ## Verify that the files are valid ##
    try:
        protein_table = pd.read_table(protein_file, low_memory=False)
    except FileNotFoundError:
        print(f"Protein file '{protein_file}' does not exist.")
        raise FileNotFoundError(f"Protein file '{protein_file}' does not exist.")
    try:
        peptide_table = pd.read_table(peptide_file, low_memory=False)
    except FileNotFoundError:
        print(f"Peptide file '{peptide_file}' does not exist.")
        raise FileNotFoundError(f"Peptide file '{peptide_file}' does not exist.")
    

    ## Create the base abundance table ##
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


    ## Rename the columns ##
    # Create a list of just the full file path names
    full_pep_run_names = peptide_cols.drop(['Precursor.Id', 'Protein.Names']).to_list()
    full_prot_run_names = protein_cols.drop(['Protein.Group', 'Protein.Names']).to_list()

    # If .raw is in the file name, strip off everything after .raw.
    if ".raw" in full_pep_run_names[0]:
        run_names = [name.split(".raw")[0] for name in full_prot_run_names]
    else:
        run_names = [os.path.splitext(name)[0] for name in full_prot_run_names]

    # Create the renaming dictionary
    pep_rename_dict = dict(zip(full_pep_run_names, run_names))
    pep_rename_dict["Precursor.Id"] = "Sequence"
    pep_rename_dict["Protein.Names"] = "Organism"

    prot_rename_dict = dict(zip(full_prot_run_names, run_names))
    prot_rename_dict["Protein.Names"] = "Organism"
    prot_rename_dict["Protein.Group"] = "Accession"

    # Rename the file columns
    pep_abundance.rename(columns=pep_rename_dict, inplace=True)
    prot_abundance.rename(columns=prot_rename_dict, inplace=True)

    # Change the organism column to be just the organism name
    pep_abundance["Organism"] = pep_abundance["Organism"].apply(lambda strg: strg.split("_")[-1])
    prot_abundance["Organism"] = prot_abundance["Organism"].apply(lambda strg: strg.split("_")[-1])

    return prot_abundance, pep_abundance


def fragpipe_22_tables(protein_file="", peptide_file="", 
                       min_unique_peptides=1, use_maxlfq=False):
    """ Process and populate pandas DataFrames with the data from the 
    FragPipe v22.0 protein and peptide files, saving only the Entry Name, 
    Peptide Sequence (for peptide), Protein ID (for protein), and the 
    file name Intensity columns. Also rename these columns.

    Parameters:
        protein_file (string): A pandas readable file with the protein data
        peptide_file (string): A pandas readable file with the peptide data
        min_unique_peptides (int): The minimum number of unique peptides
        use_maxlfq (bool): Determines whether or not to use MaxLFQ Intensity
            Default: False (uses "Intensity" columns without MaxLFQ)

    Returns:
        prot_abundance (DataFrame): pandas DataFrame with protein data
            Column names: "Organism", "Accession", and the run names
        pep_abundance (DataFrame): pandas DataFrame with peptide data
            Column names: "Organism", "Sequence", and the run names

    Raises:
        FileNotFoundError: If protein_file or peptide_file does not exist
    """
    ## Verify that the files are valid ##
    try:
        protein_table = pd.read_table(protein_file, low_memory=False)
    except FileNotFoundError:
        print(f"Protein file '{protein_file}' does not exist.")
        raise FileNotFoundError(f"Protein file '{protein_file}' does not exist.")
    try:
        peptide_table = pd.read_table(peptide_file, low_memory=False)
    except FileNotFoundError:
        print(f"Peptide file '{peptide_file}' does not exist.")
        raise FileNotFoundError(f"Peptide file '{peptide_file}' does not exist.")
    

    ## Create the base abundance table ##
    # Filters out proteins with "contam" or that are null
    peptide_table = peptide_table[~peptide_table['Protein'].str.contains("contam", na=False)]
    protein_table = protein_table[~protein_table['Protein'].str.contains("contam", na=False)]

    # Also filter out the proteins with total peptide count less than min_unique_peptides
    protein_table.query("`Combined Total Peptides` >= @min_unique_peptides", inplace=True)

    # Separate the columns of
    # "Peptide Sequence" (unique) for peptides, 
    # "Protein ID" (unique) for protein, and
    # "Entry Name" (used to filter by organism)
    peptide_cols = peptide_table.filter(regex='Peptide Sequence|Entry Name').columns.to_list()
    protein_cols = protein_table.filter(regex="Protein ID|Entry Name").columns.to_list()

    # Include the MaxLFQ columns if requested, otherwise use the the Intensity columns
    if use_maxlfq is True:
        peptide_cols += (peptide_table.filter(regex=' MaxLFQ Intensity').columns.to_list())
        protein_cols += (protein_table.filter(regex=' MaxLFQ Intensity').columns.to_list())
    else:
        # The regex uses a negative lookbehind to ignore strings with 'MaxLFQ' 
        peptide_cols += (peptide_table.filter(regex='(?<!MaxLFQ) Intensity').columns.to_list())
        protein_cols += (protein_table.filter(regex='(?<!MaxLFQ) Intensity').columns.to_list())

    # Create the abundance matrix
    pep_abundance = peptide_table.loc[:, peptide_cols]
    prot_abundance = protein_table.loc[:, protein_cols]


    # ## Rename the columns ##
    # Remove the column names that are not run names 
    peptide_cols.remove('Peptide Sequence')
    peptide_cols.remove('Entry Name')
    protein_cols.remove('Protein ID')
    protein_cols.remove('Entry Name')

    # Get just the run name (Assumes the run name has no spaces)
    run_names = [name.split(" ")[0] for name in peptide_cols]
    
    # Create the renaming dictionary
    pep_rename_dict = dict(zip(peptide_cols, run_names))
    pep_rename_dict["Peptide Sequence"] = "Sequence"
    pep_rename_dict["Entry Name"] = "Organism"
    prot_rename_dict = dict(zip(protein_cols, run_names))
    prot_rename_dict["Entry Name"] = "Organism"
    prot_rename_dict["Protein ID"] = "Accession"

    # Rename the file columns
    pep_abundance.rename(columns=pep_rename_dict, inplace=True)
    prot_abundance.rename(columns=prot_rename_dict, inplace=True)

    # Change the organism column to be just the organism name 
    # (Assumes the organism name is after an underscore "_")
    pep_abundance["Organism"] = pep_abundance["Organism"].apply(lambda strg: strg.split("_")[-1])
    prot_abundance["Organism"] = prot_abundance["Organism"].apply(lambda strg: strg.split("_")[-1])

    return prot_abundance, pep_abundance



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
    if processing_app == "diann":
        prot_abundance, pep_abundance = diann_1_9_tables(protein_file, peptide_file)
    
    elif processing_app == "fragpipe":
        prot_abundance, pep_abundance = fragpipe_22_tables(protein_file, peptide_file,
                                                           min_unique_peptides, use_maxlfq)
    
    ## Rename the DataFrames with run_IDs ##
    # Create a list of the run_IDs
    run_names = pep_abundance.columns.drop(["Sequence", "Organism"]).to_list()
    num_files = len(run_names)
    ID_generator = lambda x: str(file_id) + "-" + str(x)
    run_IDs = [ID_generator(i) for i in range(num_files)]

    # Create the table mapping run_name to run_ID
    run_metadata = pd.DataFrame({"run_name": np.array(run_names), 
                                "run_ID": np.array(run_IDs)})

    # Create a dictionary key run_name and value run_ID
    pep_rename_dict = dict(zip(run_names, run_IDs))
    prot_rename_dict = dict(zip(run_names, run_IDs))

    pep_abundance.rename(columns=pep_rename_dict, inplace=True)
    prot_abundance.rename(columns=prot_rename_dict, inplace=True)
    
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
            i += 1

        # append the information to the final_data_object for the rest
        else:
            pass
        

    return final_data_obj



def read_files(filelist = []):
    """ Read in the combined data from multiple runs. Call read_file 
    to process and populate pandas DataFrames with the data, saving 
    only the Organism, Peptide Sequence (for peptide), Protein ID 
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
    
    Raises:
        ValueError: If the processing app is not accepted
                    If the dictionary is not correct
        FileNotFoundError: If protein_file or peptide_file does not exist
    """
    # Initiate the data_objects list
    data_objects = []

    # Initiate i as a file_id counter
    i = 0

    for file_dict in filelist:
        # Give error messages if syntax is not correct
        if "protein_file" not in file_dict:
            print('Dictionary syntax incorrect. Must include "protein_file".')
            raise ValueError('Dictionary syntax incorrect. Must include "protein_file".')
        if "peptide_file" not in file_dict:
            print('Dictionary syntax incorrect. Must include "peptide_file".')
            raise ValueError('Dictionary syntax incorrect. Must include "peptide_file".')
        if "processing_app" not in file_dict:
            print('Dictionary syntax incorrect. Must include "processing_app".')
            raise ValueError('Dictionary syntax incorrect. Must include "processing_app".')
        
    
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

    return group_data(data_objects)
        



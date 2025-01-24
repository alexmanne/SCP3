import pandas as pd


def read_files(protein_file="", peptide_file=""):
    """ Read in the combined data from multiple runs. Process and populate
    pandas DataFrames with the data, saving only the Protein.Names, 
    Precursor.Id (for peptide), Protein.Group (for protein), and 
    the file name columns

    Parameters:
        protein_file (string): A pandas readable file with the protein data
        peptide_file (string): A pandas readable file with the peptide data

    Returns:
        data_obj (dict): A dictionary with the following keys/values
            "pep_abundance": A pandas DataFrame with peptide data
            "prot_abundance": A pandas DataFrame with protein data
    """
    # Create a pandas DataFrame for peptide and protein
    try:
        protein_table = pd.read_table(protein_file, low_memory=False)
    except FileNotFoundError:
        print(f"Protein file '{protein_file}' does not exist.")
        return

    try:
        peptide_table = pd.read_table(peptide_file, low_memory=False)
    except FileNotFoundError:
        print(f"Peptide file '{peptide_file}' does not exist.")
        return
    
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

    return {"pep_abundance": pep_abundance, 
            "prot_abundance": prot_abundance }
import pandas as pd


def read_files(protein_file="", peptide_file="", organisms=[None]):
    """ Read in the combined data from multiple runs. Process and populate
    pandas DataFrames with the data. Also separate according to organism by 
    assigning organized run IDs.

    Parameters:
        protein_file (string): A pandas readable file with the protein data
        peptide_file (string): A pandas readable file with the peptide data
        organisms (list): A list of organisms that appear at the end of the Protein.Names column
            default: [None]

    Returns:
        data_obj (dict): A dictionary with the following keys/values
            "run_metadata": A pandas DataFrame connecting run IDs to actual run names
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

    # Initiate i as a file_id counter
    i = 0

    for organism in organisms:
        
        # If there is an organism provided, filter the tables to only include data with that organism
        if organism:
            peptide_table = peptide_table.loc[peptide_table["Protein.Names"].str.endswith("_" + organism)]
            protein_table = protein_table.loc[protein_table["Protein.Names"].str.endswith("_" + organism)]

        # Filters out proteins with "contam" or that are null
        peptide_table = peptide_table[~peptide_table['Protein.Group'].str.contains("contam", na=False)]
        protein_table = protein_table[~protein_table['Protein.Group'].str.contains("contam", na=False)]

        # Separate the columns of just the file names and "Precursor.Id" for peptides, 
        # and "Protein.Group" for protein
        peptide_cols = peptide_table.filter(regex='\\\\|Precursor.Id').columns
        protein_cols = protein_table.filter(regex="\\\\|Protein.Group").columns

        # Create the abundance matrix
        pep_abundance = peptide_table.loc[:, peptide_cols]
        prot_abundance = protein_table.loc[:, protein_cols]

        # Now creating the runIDs






    #             # Creates a data frame of just the file names (no .raw)
    #             run_name_list = pd.DataFrame(data={"Run Names": [os.path.splitext(x)[0] for x in file_path_cols]})

    #             # Creates a new column in the data frame that is used to identify the run by the index of the data frame
    #             run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))

    #             run_name_list.to_csv("data_obj/run_name_list_pep.tsv", sep='\t')

    #             # Rename the file columns of pep_abundance to be just the run names (no .raw) 
    #             for item in [pep_abundance,pep_abundance_MS2]:
    #                 # Yes it is reduntant to make partial_column_name_mapping, but it is to use the function generate_col...
    #                 filename_to_identifier_map = dict(zip(run_name_list["Run Names"],run_name_list["Run Identifier"]))
    #                 # filename_to_filename_map = {name: name for name in run_name_list["Run Names"]}
    #                 run_name_mapping = self.generate_column_to_name_mapping(item.columns, filename_to_identifier_map)
    #                 item.rename(columns = run_name_mapping, inplace=True)

    #             # Analyze peptide ID
    #             # Create a copy of pep_abundance to manipulate
    #             pep_ID = pep_abundance.copy()
    #             pep_abundance.drop(["Precursor.Charge"], axis=1, inplace=True)

    #             # Create a list of columns, but ignore the "Annotated Sequence" column
    #             cols = [col for col in pep_ID.columns if col != 'Annotated Sequence' and col != "Precursor.Charge"]

    #             # For all of the columns
    #             for col in cols:
    #                 # If the column is not an "object" type (if it is not a string)
    #                 if pep_ID[col].dtype != 'object':
    #                     # Replace all 0's with NaN
    #                     pep_ID[col] = pep_ID[col].replace(0, np.nan)
    #                     # Replace all numerical values (floats, ints, and exponentials) to ID
    #                     pep_ID[col] = pep_ID[col].astype(str).str.replace(r"\b\d+(\.\d*)?([eE][-+]?\d+)?\b", "ID", regex=True)

    #             # Repeat the process above with MS2
    #             pep_ID_MS2 = pep_abundance_MS2.copy()
    #             pep_abundance_MS2.drop(["Precursor.Charge"], axis=1, inplace=True)
    #             cols = [col for col in pep_ID_MS2.columns if col != 'Annotated Sequence' and col != "Precursor.Charge"]
    #             for col in cols:
    #                 if pep_ID_MS2[col].dtype != 'object': 
    #                     pep_ID_MS2[col] = pep_ID_MS2[col].replace(0, np.nan)
    #                     pep_ID_MS2[col] = pep_ID_MS2[col].astype(str).str.replace(r"\b\d+(\.\d*)?([eE][-+]?\d+)?\b", "MS2", regex=True)
                
    #             # Returns a combination of the MS2 and normal data.
    #             # There is a "MS2" when a run has data for a peptide in the MS2 
    #             # There is an "MBR" when a run has data for a peptide in just the normal data
    #             pep_ID = self.combine_diann_IDs(pep_ID, pep_ID_MS2)
    #             pep_ID.drop(["Precursor.Charge"], axis=1, inplace=True)

    #             print("finish peptide")




                
    #         ## Proteins


    #             # Get the file names from the meta table
    #             protein_path_cols = protein_table_MS2.filter(regex='\\\\|Protein.Group|Protein.Names').columns

    #             # Create data frames that are just the columns with "Protein.Group", "Protein.Names" or a filename
    #             prot_abundance = protein_table.loc[:, protein_path_cols]
    #             prot_abundance_MS2 = protein_table_MS2.loc[:, protein_path_cols]

    #             # Create a list of just the file names
    #             file_path_cols = protein_table.filter(regex='\\\\').columns

    #             # Changes the columns of prot_abundance to the path/filename (no .raw)
    #             prot_abundance.columns = [os.path.splitext(os.path.basename(x))[0] if x in file_path_cols else x for x in prot_abundance.columns]
    #             prot_abundance_MS2.columns = [os.path.splitext(os.path.basename(x))[0] if x in file_path_cols else x for x in prot_abundance_MS2.columns]

    #             # Rename Protein.Group to Accession, and removes ";----" from the column
    #             prot_abundance = prot_abundance.rename(columns={'Protein.Group': 'Accession'})
    #             prot_abundance["Accession"] =  prot_abundance["Accession"].str.replace(";.*","",regex = True)
    #             prot_abundance_MS2 = prot_abundance_MS2.rename(columns={'Protein.Group': 'Accession'})   
    #             prot_abundance_MS2["Accession"] =  prot_abundance_MS2["Accession"].str.replace(";.*","",regex = True)

    #             # Create a new data frame with just the file names and the run identifier "0-n"
    #             run_name_list = pd.DataFrame(data={"Run Names": [os.path.splitext(os.path.basename(x))[0] for x in file_path_cols]})
    #             run_name_list['Run Identifier'] = run_name_list.index.to_series().apply(lambda x: str(file_id) + "-" + str(x))

    #             run_name_list.to_csv("data_obj/run_name_list_prot.tsv", sep='\t')

    #             for item in [prot_abundance, prot_abundance_MS2]:
    #                 # Generate a new column name mapping using the function
    #                 filename_to_identifier_map = dict(zip(run_name_list["Run Names"], run_name_list["Run Identifier"]))
    #                 # filename_to_filename_map = {name: name for name in run_name_list["Run Names"]}
    #                 fileid_mapping = self.generate_column_to_name_mapping(item.columns, filename_to_identifier_map)
    #                 item.rename(columns=fileid_mapping, inplace=True)


    #         # Analyze Protein ID
    #             # Create a copy of prot_abundance to manipulate
    #             prot_ID = prot_abundance.copy()
    #             prot_abundance.drop(["Protein.Names"], axis=1, inplace=True)

    #             # Create a list of columns, but ignore the "Accession" column
    #             cols = [col for col in prot_ID.columns if col != 'Accession' and col != "Protein.Names"]

    #             # For all of the columns
    #             for col in cols:
    #                 # If the column is not an "object" type (if it is not a string)
    #                 if prot_ID[col].dtype != 'object':
    #                     # Replace all 0's with NaN
    #                     prot_ID[col] = prot_ID[col].replace(0, np.nan)
    #                     # Replace all numerical values (floats, ints, and exponentials) to ID
    #                     prot_ID[col] = prot_ID[col].astype(str).str.replace(r"\b\d+(\.\d*)?([eE][-+]?\d+)?\b", "ID", regex=True)

    #             # Repeat the the process above with MS2
    #             prot_ID_MS2 = prot_abundance_MS2.copy()
    #             prot_abundance_MS2.drop(["Protein.Names"], axis=1, inplace=True)
    #             cols = [col for col in prot_ID_MS2.columns if col != 'Accession' and col != "Protein.Names"]
    #             for col in cols:
    #                 if prot_ID_MS2[col].dtype != 'object':
    #                     prot_ID_MS2[col] = prot_ID_MS2[col].replace(0, np.nan)
    #                     prot_ID_MS2[col] = prot_ID_MS2[col].astype(str).str.replace(r"\b\d+(\.\d*)?([eE][-+]?\d+)?\b", "MS2", regex=True)


    #             # Returns a combination of the MS2 and normal data.
    #             # There is a "MS2" when a run has data for a peptide in the MS2 
    #             # There is an "MBR" when a run has data for a peptide in just the normal data
    #             prot_ID = self.combine_diann_IDs(prot_ID, prot_ID_MS2)   
    #             prot_ID.drop(["Protein.Names"], axis=1, inplace=True)

    #             print("finish protein")


    #     # Creates a column in run_name_list for the processing app 
    #     run_name_list["Processing App"] = process_app

    #     run_name_list.to_csv("data_obj/run_name_list.tsv", sep='\t')
    #     return_matrix = {'run_metadata': run_name_list}  
    #     if self.ignore_peptides is False:
    #         peptide_ID_summary = self.sumIDs(pep_ID)  
    #         return_matrix['peptide_other_info']=  pep_other_info
    #         return_matrix['peptide_abundance']=  pep_abundance
    #         return_matrix['peptide_ID_matrix']=  pep_ID
    #         return_matrix['peptide_ID_Summary']=  peptide_ID_summary

    #     if self.ignore_proteins is False:
    #         protein_ID_summary = self.sumIDs(prot_ID)
    #         return_matrix['protein_other_info']= prot_other_info
    #         return_matrix['protein_abundance']=  prot_abundance
    #         return_matrix['protein_ID_matrix']=  prot_ID
    #         return_matrix['protein_ID_Summary']=  protein_ID_summary

    #     return return_matrix




    # # Initiate a list of data objects
    # data_objects = []

    # # Initiate i as a file_id counter
    # i = 0

    # # Loop through the organisms in organisms for the sake of creating the file_id
    # for organism in organisms:
    #     data_object = read_file(protein_file, peptide_file, 
    #                             organism, file_id=i)
    #     data_objects.append(data_object)


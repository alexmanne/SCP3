
def analyze_data(data_obj, settings):
    """ Compute the statistical analysis of the objects in data_obj 
    based on the settings. Executes the following analyses.

    - Filter by Condition
    - Filter by Missing Value
    - Filter by Unique Peptide (Happens in read_files)
    - Normalize to Median

    Parameters:
        data_obj (list): A list of dictionaries that contain the data objects 
            to combine. Each data object has the following keys/values     
                "run_metadata": A pandas DataFrame mapping run names to ids    
                "pep_abundance": A pandas DataFrame with peptide data     
                "prot_abundance": A pandas DataFrame with protein data    
        settings (dict): A dictionary containing the settings from the      
            settings file

    Returns:

    
    """


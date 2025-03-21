

def filter_by_missing_value(abundance_data, missing_value_threshold=33, 
                            is_protein=True):
    """ Filter our proteins that have a higher percent of missing values 
    than the threshold

    Parameters:
        abundance_data (DataFrame): A Pandas DataFrame containing the
            abundance data to filter
        missing_value_threshold (int): An integer outlining the threshold.   
        is_protein (bool): Tells if this is protein or peptide data

    Returns:
        abundance_data (DataFrame): The filtered data
    
    """
    if is_protein:
        identifying_col = "Accession"
    else:
        identifying_col = "Sequence"

    # Get the run columns
    data_columns = abundance_data.columns.drop([identifying_col, "Protein Name"]).to_list()
    num_columns = len(data_columns)

    # Create a column counting the missing values in each row
    abundance_data["missing_values"] = abundance_data.isna().sum(axis=1)

    # Turn that count into a percentage of the runs
    abundance_data["missing_value_rate"] = abundance_data["missing_values"] / num_columns * 100
    abundance_data = abundance_data.query(
            "missing_value_rate <= @missing_value_threshold")

    return abundance_data.drop(columns=["missing_values", "missing_value_rate"])
    

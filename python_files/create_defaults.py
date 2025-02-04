

def generate_input_files_settings(n):
    """ Writes a yaml file with empty values for n files 

    Parameters:
        n (int): number of files
    """
    comment = ("# This file creates a list of the files to import\n"
               "\n"
               "# Be sure to include the following key words\n"
               "# - peptide_file: 'path_to_file'\n"
               "#   protein_file: 'path_to_file'\n"
               "#   processing_app: diann\n"
               "\n"
               "# Use a '-' to start a new file group\n\n\n")
    
    data = ("- peptide_file: input/PUT FILE NAME HERE\n"
            "  protein_file: input/PUT FILE NAME HERE\n"
            "  processing_app: diann    # Options are diann or fragpipe\n")

    with open("default_file_settings.yaml", "w") as file:
        file.write(comment)

        for i in range(n):
            file.write(f"# File {i}\n")
            file.write(data + "\n")


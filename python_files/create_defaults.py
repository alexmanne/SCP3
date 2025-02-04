

def generate_input_files_template(n):
    """ Writes a yaml file with empty values for n files.

    Parameters:
        n (int): number of files

    Raises:
        FileExistsError: If the name "default_file_settings.yaml" already exists
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

    # Use a Try except block to signify if the file already exists
    try:
        with open("default_file_template.yaml", "x") as file:
            file.write(comment)

            for i in range(n):
                file.write(f"# File {i}\n")
                file.write(data + "\n")
    except FileExistsError:
        print("A file with the name 'default_file_template.yaml' already exists.\n" 
              "Change the name of your file before generating the template.")
        raise FileExistsError("Will not overwrite 'default_file_template.yaml'.")


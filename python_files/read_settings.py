import yaml


def read_settings(settings_file):
    """ Read in the settings yaml file to save them for the
    rest of the run.

    Parameters:
        settings_file (string): A string name for the settings yaml file

    Return:
        settings_dict (dict): A dictionary with all of the saved settings

    Raises:
        FileNotFoundError: If the yaml_file does not exist
    """
    ## Verify that the file is valid ##
    try:
        with open(settings_file, "r") as f:
            settings_dict = yaml.safe_load(f)
    except FileNotFoundError:
        print(f"File '{settings_file}' does not exist.")
        raise FileNotFoundError(f"File '{settings_file}' does not exist.")
    
    return settings_dict


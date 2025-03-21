import sys
sys.path.insert(0, '/Users/genemannewitz/Documents/KellyLab/SCP3')

from python_files.filter_data import filter_by_missing_value
import pandas as pd

def test_filter_by_missing_value():
    filename = "/Users/genemannewitz/Documents/KellyLab/SCP3/testing/small_dia_report.pg_matrix.csv"
    small_abundance = pd.read_csv(filename)

    # Test for missing value threshold = 33
    filtered_abundance = filter_by_missing_value(small_abundance, 33)
    proteins_left = filtered_abundance.shape[0]
    assert(6, proteins_left), f"Number of proteins left is {proteins_left}"

    # Test for missing value threshold = 40
    filtered_abundance = filter_by_missing_value(small_abundance, 40)
    proteins_left = filtered_abundance.shape[0]
    assert(7, proteins_left), f"Number of proteins left is {proteins_left}"

    # Test for missing value threshold = 66
    filtered_abundance = filter_by_missing_value(small_abundance, 66)
    proteins_left = filtered_abundance.shape[0]
    assert(8, proteins_left), f"Number of proteins left is {proteins_left}"

    # Test for missing value threshold = 100
    filtered_abundance = filter_by_missing_value(small_abundance, 100)
    proteins_left = filtered_abundance.shape[0]
    assert(18, proteins_left), f"Should not have filtered any out"

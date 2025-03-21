
def generate_input_files_template(n):
    """ Writes a yaml file with empty values for n files.

    Parameters:
        n (int): number of files

    Raises:
        FileExistsError: If the name "input_file_template.yaml" already exists
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
            "  processing_app: diann    # Options: diann, fragpipe, spectronaut\n")

    # Use a Try except block to signify if the file already exists
    try:
        with open("input_file_template.yaml", "x") as file:
            file.write(comment)

            for i in range(n):
                file.write(f"# File {i}\n")
                file.write(data + "\n")
    except FileExistsError:
        print("A file with the name 'input_file_template.yaml' already exists.\n" 
              "Change the name of your file before generating the template.")
        raise FileExistsError("Will not overwrite 'input_file_template.yaml'.")



def generate_settings_template():
    """ Writes a yaml file with default settings.

    Raises:
        FileExistsError: If the name "settings_template.yaml" already exists
    """

    default = """# Proteomics Pipeline Settings File

# Experimental Design
experimental_design:
  conditions:
    - name: Condition_A
      replicates: [Sample_1, Sample_2, Sample_3]
    - name: Condition_B
      replicates: [Sample_4, Sample_5, Sample_6]
  control_condition: Condition_A

  # Randomization and Blocking
  randomization:
    enable: true  # Enable randomization
    seed: 42  # Random seed for reproducibility
    method: complete_randomization  # Options: complete_randomization, stratified, etc.

  blocking:
    enable: true  # Enable blocking
    block_by: batch  # Variable to block by (e.g., batch, day, etc.)
    blocks:
      - name: Batch_1
        samples: [Sample_1, Sample_4]
      - name: Batch_2
        samples: [Sample_2, Sample_5]
      - name: Batch_3
        samples: [Sample_3, Sample_6]

# Data Filters
filters:
  min_peptides: 2      # fragpipe only
  min_confidence_score: 0.95
  min_peptide_length: 7
  max_missing_values: 0.5
  use_maxlfq: False    # fragpipe only

# Normalization Options
normalization:
  method: median_normalization  # Options: median_normalization, quantile_normalization, none
  reference_condition: Condition_A

# Log Transformation
log_transformation:
  apply: true
  base: 2  # Options: 2, 10

# Statistical Analysis Settings
statistical_analysis:
  # p-value settings
  p_value:
    use_p_value: true  # Set to true to use p-values
    p_value_threshold: 0.05  # Threshold for significance (p < 0.05)
  
  # q-value (FDR-adjusted p-value) settings
  q_value:
    use_q_value: false  # Set to true to use q-values
    q_value_threshold: 0.05  # Threshold for significance (q < 0.05)
    false_discovery_rate: 0.01  # FDR threshold for multiple testing correction

  # Confidence interval
  confidence_interval: 0.95  # Confidence level for statistical tests

## NEW PAGE ##
# Plotting Options
plotting:
  # General settings for all plots
  general:
    theme: classic  # Options: classic, minimal, dark, etc.
    color_palette: viridis  # Options: viridis, plasma, magma, inferno, cividis, etc.
    font_size: 12  # Base font size for all plots
    output_directory: "./plots"  # Directory to save plots
    file_format: png  # Options: png, pdf, svg, jpeg
    dpi: 300  # Resolution for saved plots

  # ID Plot (Identification Plot)
  id_plot:
    enable: true
    x_axis: sample_id  # X-axis variable (e.g., sample IDs)
    y_axis: protein_count  # Y-axis variable (e.g., number of proteins identified)
    title: "Protein Identification Plot"
    x_label: "Sample ID"
    y_label: "Number of Proteins Identified"
    color_by: condition  # Color points by condition
    show_labels: true  # Show sample labels on the plot

  # CV Violin Plot (Coefficient of Variation)
  cv_violin_plot:
    enable: true
    group_by: condition  # Group violins by condition
    title: "Coefficient of Variation (CV) Distribution"
    x_label: "Condition"
    y_label: "CV (%)"
    show_points: true  # Overlay individual data points
    point_size: 2  # Size of overlaid points

  # Venn Diagram
  venn_diagram:
    enable: true
    groups:  # List of groups to compare
      - Condition_A
      - Condition_B
      - Condition_C
    title: "Venn Diagram of Protein Overlap"
    colors:  # Custom colors for each group
      Condition_A: "#1f77b4"
      Condition_B: "#ff7f0e"
      Condition_C: "#2ca02c"

  # Volcano Plot
  volcano_plot:
    enable: true
    x_axis: log2_fold_change  # X-axis variable
    y_axis: -log10_p_value  # Y-axis variable
    title: "Volcano Plot"
    x_label: "Log2 Fold Change"
    y_label: "-Log10 p-value"
    significance_threshold: 0.05  # Threshold for significance (p-value or q-value)
    log2_fc_threshold: 1  # Log2 fold change threshold
    color_significant: "#d62728"  # Color for significant points
    color_nonsignificant: "#2ca02c"  # Color for non-significant points
    show_labels: true  # Label significant points
    label_threshold: 0.01  # Threshold for labeling points

  # Heatmap
  heatmap:
    enable: true
    data: protein_abundance  # Data to plot (e.g., protein abundance matrix)
    title: "Protein Abundance Heatmap"
    x_label: "Samples"
    y_label: "Proteins"
    clustering_method: complete  # Options: complete, average, single
    color_scheme: viridis  # Options: viridis, plasma, magma, inferno, cividis
    show_row_dendrogram: true  # Show dendrogram for rows (proteins)
    show_col_dendrogram: true  # Show dendrogram for columns (samples)

  # PCA (Principal Component Analysis)
  pca:
    enable: true
    data: protein_abundance  # Data to plot (e.g., protein abundance matrix)
    title: "PCA Plot"
    color_by: condition  # Color points by condition
    shape_by: replicate  # Shape points by replicate
    show_labels: true  # Show sample labels
    show_ellipses: true  # Show confidence ellipses for each group
    ellipse_confidence: 0.95  # Confidence level for ellipses

  # GSEA (Gene Set Enrichment Analysis)
  gsea:
    enable: true
    gene_sets:  # List of gene sets to analyze
      - GO_Biological_Process
      - KEGG_Pathways
    title: "GSEA Enrichment Plot"
    x_label: "Rank in Ordered Dataset"
    y_label: "Enrichment Score"
    color_scheme: plasma  # Color scheme for enrichment plots
    show_leading_edge: true  # Highlight leading edge genes

  # Missing Value Plot
  missing_value_plot:
    enable: true
    data: protein_abundance  # Data to plot (e.g., protein abundance matrix)
    title: "Missing Value Distribution"
    x_label: "Samples"
    y_label: "Proteins"
    color_by: condition  # Color bars by condition
    show_percentage: true  # Show percentage of missing values

# Output Options
output:
  directory: "./results"
  file_format: csv  # Options: csv, tsv, xlsx
  include_raw_data: true
  include_normalized_data: true
  include_statistical_results: true
"""
    
    # Use a Try except block to signify if the file already exists
    try:
        with open("settings_template.yaml", "x") as file:
            file.write(default)

    except FileExistsError:
        print("A file with the name 'settings_template.yaml' already exists.\n" 
              "Change the name of your file before generating the template.")
        raise FileExistsError("Will not overwrite 'settings_template.yaml'.")



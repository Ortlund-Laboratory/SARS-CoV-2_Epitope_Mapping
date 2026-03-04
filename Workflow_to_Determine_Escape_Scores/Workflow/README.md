# Workflow to Calculate Escape Scores

SARS-CoV-2 nucleocapsid variant escape scores were investigated, with a reference and a low-binding population collected for sequencing per antibody. These sequences can be processed to give escape scores for each variant. We use Ab339 from CorDx as an example here.

## Input Files Required

[SnakeFile](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/Snakefile)<br>
Gives overall instructions for the `snakemake` workflow.<br>
[config.yaml](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/config.yaml)<br>
Configuration script controlling variables used by Jupyter notebooks.<br>
[build_variants.ipynb](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/build_variants.ipynb)<br>
Builds a barcode variant table based on the data from the processed PacBio CCSs.<br>
[R2_to_R1.py](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/R2_to_R1.py)<br>
Converts barcodes located at the R2 end to the R1 end by taking the reverse complement. This allows the barcodes to be read and parsed correctly by the [illuminabarcodeparser](https://jbloomlab.github.io/dms_variants/dms_variants.illuminabarcodeparser.html#dms_variants.illuminabarcodeparser.IlluminaBarcodeParser) algorithm.<br>
[count_variants.ipynb](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/count_variants.ipynb)<br>
Counts the number of times a barcode (and by extension a variant) appears in each Illumina barcode sequencing sample.<br>
[scripts/run_nb.py](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/scripts/run_nb.py)<br>
Runs Jupyter notebooks and creates Markdown output.<br>
[data/feature_parse_specs.yaml](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/data/feature_parse_specs.yaml)<br>
Script for controlling the sequence parsing strategy.<br>
[data/PacBio_amplicons.gb](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/data/PacBio_amplicons.gb)<br>
GeneBank data file describing sequence features.<br>
[data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/data/barcode_runs.csv)<br>
List of Illumina barcode samples to be analyzed by the snakemake workflow.<br>
[data/processed_ccs.csv.gz](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/data/processed_ccs.csv.gz)<br>
Processed PacBio CCSs, generated from our [PacBio](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/PacBio) routine. Ensure the library is consistent with those used for the assay.<br>
[data/wildtype_sequence.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/data/wildtype_sequence.fasta)<br>
SARS-CoV-2 nucleocapsid wildtype sequence (Wuhan).<br>

### Sequencing Data

The workflow operates on Illumina barcode sequencing data in fastq.gz format and these files are kept compressed throughout. File location and name should match the listings given in [data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/data/barcode_runs.csv). These files are too large to be contained in GitHub, and will be uploaded to the NCBI along with nucleocapsid variant sequences following FACS analysis with all the other antibodies we investigated. For reviewers, who may wush for access in the meantime, the sequences are available upon request (adkeith@emory.edu & eortlun@emory.edu).

## Workflow

Use the `snakemake` environment:

`conda activate snakemake`

Run `snakemake` using specified number of cores:

`snakemake -j 6`

## Key Output

[results/counts/barcode_fates.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/blob/main/Workflow_to_Determine_Escape_Scores/Workflow/results/counts/barcode_fates.csv)<br>
Tally of barcodes classified and filtered according to quality.<br>
[results/counts/variant_counts.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/blob/main/Workflow_to_Determine_Escape_Scores/Workflow/results/counts/variant_counts.csv)<br>
Tally of individual barcode counts for each sample.<br>

## Stability Score Generation & Data Visualization

Go to [scores_and_visualization](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/tree/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization) for higher-/lower-stability score calculation and heatmap generation.

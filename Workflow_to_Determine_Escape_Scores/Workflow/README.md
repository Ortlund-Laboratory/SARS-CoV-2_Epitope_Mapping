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

The workflow operates on Illumina barcode sequencing data in fastq.gz format and these files are kept compressed throughout. File location and name should match the listings given in [data/barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/data/barcode_runs.csv). These files are too large to be contained in GitHub, and are available from the NCBI, along with nucleocapsid variant sequences following FACS analysis with all the other antibodies we investigated.

NCBI links are as follows:
BioProject:
BioSample:

Please consult [Log_of_files_for_NCBI.xlsx]() to determine the appropriate reference files to be used when replicating escape score results. We always collected a reference population on the same day we ran FACS experiments to ensure consistency between cells. Please note that a reference file for the R001_Quidel antibody is missing, and sequencing datasets for MAD4904_Roche, mAb22-048-C_Roche, mAb-OTI-1_OraSure, mAb-OTI-2_OraSure, mAb4_AzureBiotech and mAb5_AzureBiotech are no longer available, as detailed in the Excel file. Though these sequencing files are not available, we have retained the escape score files (see [here](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Individual_Antibody_Escape_Profiles) and so they can be included in downstream analyses. 


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

## Escape Score Generation

Go to [escape_score_calc](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/escape_score_calc) for escape score calculation.

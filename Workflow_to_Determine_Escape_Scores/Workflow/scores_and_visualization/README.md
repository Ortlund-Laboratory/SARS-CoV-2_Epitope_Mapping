# File Prep

[variant_counts.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/results/counts/variant_counts.csv) from results/counts needs to be manipulated.

By referring to [barcode_runs.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/data/barcode_runs.csv), we can see that two experiments were run: a reference sample (exp1-none-0-reference) and the escape population (exp2-Ab339-2000-escape).

For compatibility with our R scripts, we require the information in [variant_counts.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/results/counts/variant_counts.csv) to be split into three-column (barcode,mutation,count) files. This is very simple and you can write a script to do this if you like.

First, remove all variants where the number of amino acid mutations is 0 (i.e. is still wildtype) or greater than 1 (multiple point mutations).

```
awk -F',' '!($10!=1)' variant_counts.csv > tmp.csv && mv tmp.csv variant_counts.csv
```
Next step is to separate into reference, mycpos and mycneg files.
```
sed -n '/exp1-none-0-reference/p' variant_counts.csv > ref_variant_counts.txt
sed -n '/exp2-Ab339-2000-escape/p' variant_counts.csv > esc_variant_counts.txt
```
Then reformat each file so that they just contain tab-separated columns for barcode, mutation and count.
```
awk '{print $4, $8, $5}' FS="," OFS="\t" ref_variant_counts.txt > tmp.txt && mv tmp.txt ref_variant_counts.txt
awk '{print $4, $8, $5}' FS="," OFS="\t" esc_variant_counts.txt > tmp.txt && mv tmp.txt esc_variant_counts.txt
```
Finally, add tab-separated headers (barcode  mutation  count) for each of these files.

Now that we have the input files formatted correctly, we are ready to calculate escape scores and to visualize the data.

# Protocol

## Input Files Required

[Mar2022_mycneg_myBarcodeMapping.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/Mar2022_mycneg_myBarcodeMapping.R)<br>
R script to calculate scores from reference and mycneg counts, generate heatmaps and produce files to map scores onto structures.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[mycneg_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/mycneg_variant_counts.txt)<br>
Barcode counts for the mycneg sample.<br>
[N_Wuhan.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/N_Wuhan.fasta)<br>
Amino acid sequence for N Wuhan. This is required to complete the heatmap.

## Workflow

```
rstudio Mar2022_mycneg_myBarcodeMapping.R
```

## Key Output

[MycNeg_scores.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/output/MycNeg_scores.csv)<br>
Log of each mutation and its associated score.<br>
[MycNeg_per_site_score.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/output/MycNeg_per_site_score.csv)<br>
Log of each site and its associated averaged score.<br>
[MycNeg_Fraction_heatmap01.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/output/MycNeg_Fraction_heatmap01.png)<br>
Heatmap of scores for all single-point variants from the MycNeg sample (part 1).<br>
[MycNeg_Fraction_heatmap02.png](https://github.com/Ortlund-Laboratory/SARS-CoV-2-Structure/blob/main/Raw%20Deep%20Mutational%20Scanning%20(DMS)%20Data/Workflow/scores_and_visualization/mycneg/full_NP/output/MycNeg_Fraction_heatmap02.png)<br>
Heatmap of scores for all single-point variants from the MycNeg sample (part 2).<br>


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

Now that we have the input files formatted correctly, we are ready to calculate escape scores.

# Protocol

## Input Files Required

[BarcodeMapping.R](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/blob/main/Workflow_to_Determine_Escape_Scores/Workflow/escape_score_calc/BarcodeMapping.R)<br>
R script to calculate escape scores from reference and escape counts.<br>
[ref_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/escape_score_calc/ref_variant_counts.txt)<br>
Barcode counts for the reference sample.<br>
[esc_variant_counts.txt](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/escape_score_calc/esc_variant_counts.txt)<br>
Barcode counts for the escape sample.<br>
[N_Wuhan.fasta](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/escape_score_calc/N_Wuhan.fasta)<br>
Amino acid sequence for N Wuhan.

## Workflow

```
rstudio BarcodeMapping.R
```

## Key Output

[Ab339_escape_fractions.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/escape_score_calc/Ab339_escape_fractions.csv)<br>
Log of each mutation and its associated score.<br>
[Ab339_average_escape.csv](https://github.com/Ortlund-Laboratory/SARS-CoV-2_Epitope_Mapping/tree/main/Workflow_to_Determine_Escape_Scores/Workflow/escape_score_calc/Ab339_average_escape.csv)<br>
Log of each site and its associated averaged score.<br>


The output here is not z-normalized, but can be for further analysis.


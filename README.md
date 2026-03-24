# Bacterial Motility Analysis

Welcome!

This script was developed as part of my PhD project to analyse **bacterial motility on low-agar-content plates**. It is designed to work with a separate R file containing the experimental data matrix and provides both **publication-ready visualization** and **basic statistical comparisons** of motility diameter measurements.

---

## Overview

The script analyses motility assays for two bacterial strains measured under multiple treatment conditions and in association with two spider species:

- **PRD** – *Pardosa lugubris*
- **PTSD** – *Parasteatoda tepidariorum*

The workflow includes:

- conversion of the input matrix into long format,
- generation of faceted boxplots with jittered individual observations,
- comparison of selected treatments against the control group,
- multiple-testing correction,
- annotation of statistical significance directly on the plot,
- export of the final figure as a high-resolution TIFF file.

---

## Experimental Design

The script was written for a dataset with the following structure:

- **12 replicates**
- **4 treatments**
- **2 spider species**

The measured variable is:

- **motility diameter (cm)**

The default treatment groups are:

- `CONTROL`
- `ANTIBIO`
- `EGGS`
- `SILK`

The default bacterial strains recognized by the script are:

- *Escherichia coli* (NCTC 9001)
- *Pseudomonas aeruginosa* (NCTC 10662)

---

## Input Requirements

To run this script, you need a separate R file containing the data matrix, for example:

`source("/Users/Mateusz/OneDrive/Desktop/motility_data_matrices.R")`

This file must create an object named:

`M`
(If the object M does not exist after sourcing the file, the script will stop with an error.)

## Expected data structure

The script assumes that:

- rows represent biological replicates,
- columns represent combinations of spider species and treatment groups,
- row names encode bacterial strain identity.

The matrix should follow a structure compatible with:

- 12 replicates × 4 treatments × 2 spider species

## Statistics
By default, the script uses:

- Wilcoxon rank-sum test
- Holm correction for multiple comparisons

Comparisons are performed separately within each strain × spider species combination.

## Output

The script generates:

- a faceted boxplot with jittered points,
- a summary table of statistical results,
- a final annotated plot with significance stars.

## Notes

This script was written to be adaptable for future motility experiments with a similar design.
Some internal comments may still contain Polish.

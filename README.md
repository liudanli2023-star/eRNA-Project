# eRNA-project

This repository contains scripts used for enhancer RNA (eRNA) analysis during mouse brain development. The workflow supports construction of enhancer-gene candidate regulatory pairs, annotation of developmental stages, integration with marker genes, screening of mouse-human conserved enhancer regions, and annotation with Hi-C loop evidence.

## Workflow

1. `scripts/step1_matrix.R`
   - Calculates Pearson correlation coefficients and corresponding P values between enhancer TPM profiles and gene TPM profiles across developmental stages.
   - Outputs correlation and P-value matrices in RData format.

2. `scripts/step2.pl`
   - Generates enhancer-gene candidate pairs based on genomic distance.
   - Uses enhancer coordinates and gene annotation to identify genes located within 1 Mb of each enhancer center.

3. `scripts/step3_1.R`
   - Performs FDR correction and filters significant enhancer-gene pairs.
   - Calculates a significance score based on correlation and corrected P value.

4. `scripts/step3_2.R`
   - Assigns the developmental stage with maximum enhancer expression to each significant enhancer-gene pair.
   - Also summarizes cases where one enhancer is associated with multiple candidate genes.

5. `scripts/step3_3.R`
   - Adds enhancer TPM and gene TPM values at the corresponding developmental stage.
   - Produces tab-delimited and CSV outputs for downstream analysis.

6. `scripts/step4_markergene.R`
   - Screens enhancer-gene regulatory pairs associated with predefined mouse brain developmental marker genes.
   - Adds cell-type annotation and applies expression/length filtering.

7. `scripts/Step4_orthologs.R`
   - Screens mouse enhancer-gene pairs whose enhancer regions have cross-species orthologous coordinates in the human genome.

8. `scripts/step4_ASDgene.R`
   - Intersects conserved enhancer-gene pairs with autism spectrum disorder (ASD)-related genes.

9. `scripts/Step5_addLoop.R`
   - Adds Hi-C loop support information to marker gene-associated enhancer-gene pairs.
   - Reports full enhancer-gene matches and enhancer-matched/gene-different relationships.

10. `scripts/Step5_ASD_addLoop.R`
    - Adds Hi-C loop support information to ASD gene-associated enhancer-gene pairs.

## Main Inputs

- Enhancer TPM matrix.
- Gene TPM matrix.
- Enhancer genomic coordinates.
- Gene annotation file in GFF3/GTF format.
- Mouse-human enhancer orthologous coordinate file.
- Marker gene list or marker gene-cell type table.
- Hi-C promoter-enhancer loop files.
- ASD risk gene list.

## Notes

- Several scripts contain absolute paths from the original analysis environment. Before running the workflow on another machine, update `setwd()` and input file paths according to the local directory structure.
- Large intermediate matrices are saved in RData format to reduce repeated calculation.
- Thresholds such as correlation coefficient, P value, enhancer TPM, and enhancer length can be adjusted according to the study design.

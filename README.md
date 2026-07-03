# EEG-psychotropic-medications
Code used for automatic EEG data preprocessing and calculating main effects described in publication  "Brainwaves Under Medication: Revealing Class-Specific Neural Signatures of Psychotropic Medication from 24,000 EEGs"

The matlab file "Brainwaves_Under_Medication_code_commented" calculates effect sizes and p values for PCA components or EEG features for different medicine groups (for two comparisons, one-vs-rest, Ovr, and with drug-naive subjects, DN). P-values are corrected for multiple comparisons using bonf_holm helper function.
The code uses helper functions balance_groups_meds to create comparison groups with similat age, sex and disorder distributions.

The python notebook "Brainwaves_Under_Medication_Visualization_and_PCA_analysis" creates graphs for visualizationg these results, including exemplificatory PCA components, features contributing the most to significant PCA components for the drug classes, and mean feature values comparisons between the groups.

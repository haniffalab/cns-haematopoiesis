Metadata overview

This directory contains curated metadata tables describing samples used in the Bone Atlas and in cross-organ developmental comparison datasets. These metadata files provide harmonised information on sample origin, technology, age, and data access, and are used throughout downstream analyses.

Core combined metadata
sample_metadata/COMBINED_BONE_CROSS_ORGAN_METADATA_FINAL.csv

This file contains harmonised metadata for:

Bone Atlas samples

Developmental organ datasets used for cross-organ comparison

Each row corresponds to a single 10X library.

Key columns

Sample_ID
10X sequencing lane or library identifier

Organ
Broad organ or tissue category (e.g. bone, brain, liver)

Anatomical_Site
More specific anatomical sampling location

10X_Chemistry
10X Genomics chemistry version used

Method
Library preparation method (5′, 3′, or multiome)

Technology
Cellular or nuclear RNA sequencing

Age
Developmental age of the sample

Multiplexed
Whether samples were multiplexed (Yes / No)

Study
Name of the publication or study associated with the dataset

Accession_ID
Raw data accession or storage location (e.g. GEO, EGA, ArrayExpress)

FACS_Sorting
Whether the sample was FACS-sorted prior to sequencing (Yes / No)

Organ-specific metadata files
CROSS_ORGAN_METADATA.csv

Metadata for developmental organ datasets only, excluding Bone Atlas samples.

Linnarson_EGAD00001006049.csv

Metadata for the brain dataset from
Braun et al., 2022
(EGA accession: EGAD00001006049)

BONE_ATLAS_METADATA.csv

Metadata for Bone Atlas samples only, covering bone-derived organs and anatomical sites.

ArrayExpress / BioStudies metadata
cross_organ_comparison_EMAT_metadata/

This directory contains metadata files for developmental organ datasets used in cross-organ comparisons. These were obtained from BioStudies / ArrayExpress:

Study accession: E-MTAB-14385

Link: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14385

Files correspond to original E-MTAB metadata tables.

bone_atlas_EMAT_metadata/

This directory contains metadata files for Bone Atlas datasets obtained from BioStudies / ArrayExpress:

Study accession: E-MTAB-14385

Link: https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14385

These files represent the original submitted metadata for Bone Atlas samples.
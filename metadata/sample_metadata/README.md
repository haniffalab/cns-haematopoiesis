# Metadata directory

This directory contains metadata tables used to describe samples included in the Bone Atlas and in cross-organ developmental comparison datasets.  
Files are organised by scope (combined, bone-specific, cross-organ, and legacy/reference metadata).

---

## `sample_metadata/`

Metadata tables that are actively used for analysis and harmonisation.

### `COMBINED_BONE_CROSS_ORGAN_METADATA_FINAL.csv`

Harmonised metadata for all samples included in the Bone Atlas and in developmental datasets used for cross-organ comparison.

- Each row corresponds to a single 10X library.
- This is the **primary metadata file** used throughout the project.
- Includes information on sample origin, anatomical site, age, technology, and raw data accessions.

Key fields include:
- Sample ID / 10X library ID  
- Organ and anatomical site  
- Developmental age  
- Sequencing chemistry and method  
- Cellular vs nuclear RNA  
- Associated study and accession ID  
- Multiplexing and FACS sorting status  

---

## `CROSS_ORGAN_METADATA.csv`

Metadata for **developmental organ datasets only**, excluding Bone Atlas samples.

- Used for analyses restricted to non-bone organs.
- Subset of the combined metadata with bone-specific samples removed.

---

## `BONE_ATLAS_METADATA.csv`

Metadata for **Bone Atlas samples only**.

- Includes all bone-derived organs and anatomical sites.
- Used for bone-restricted analyses and figure generation.

---

## Study-specific metadata files

### `Linnarson_EGAD00001006049.csv`

Metadata for the brain dataset from **Braun et al., 2022**.

- EGA accession: `EGAD00001006049`
- Contains study-specific annotations provided with the original publication.

---

## `cross_organ_comparison_EMAT_metadata/`

Original metadata files for **developmental organ datasets** obtained from BioStudies / ArrayExpress.

e.g.; 
- https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14385

These files reflect the metadata as deposited in public repositories and are retained for reference and traceability.

---

## `bone_atlas_EMAT_metadata/`

Original metadata files for **Bone Atlas datasets** obtained from BioStudies / ArrayExpress.

--

<!--- This markdown file was designed to roughly follow the Penn LINC Neuroinformatics template: https://pennlinc.github.io/docs/Contributing/ProjectTemplate/ --->

# Predicting future cognitive impairment in preclinical Alzheimer’s disease using amyloid PET and MRI: a multisite machine learning study

![Graphical abstract](figures/figure1_overview.png)

# Project Description

This repository contains accompanying code for the manuscript "[Predicting future cognitive impairment in preclinical Alzheimer’s disease using amyloid PET and MRI: a multisite machine learning study](https://www.medrxiv.org/content/10.1101/2025.10.15.25337507v2)".

## Corresponding Authors

- [Braden Yang](mailto:b.y.yang@wustl.edu)
- [Aris Sotiras](mailto:aristeidis.sotiras@pennmedicine.upenn.edu)

## Coauthors

- Tom Earnest
- Murat Bilgel
- Marilyn S Albert
- Sterling C Johnson
- Christos Davatzikos
- Guray Erus
- Colin L Masters
- Susan M Resnick
- Michael I Miller
- Arnold Bakker
- John C Morris
- Tammie Benzinger
- Brian Gordon

## Datasets

- [Anti-Amyloid Treatment in Asymptomatic Alzheimer’s (A4) Study](https://www.a4studydata.org/)
- [Alzheimer's Disease Neuroimaging Initiative (ADNI)](https://adni.loni.usc.edu/)
- [Harvard Aging Brain Study (HABS)](https://habs.mgh.harvard.edu/researchers/)
- [Mayo Clinic Study of Aging (MCSA)](https://ida.loni.usc.edu/collaboration/access/appLicense.jsp)
- [Preclinical Alzheimer's Disease Consortium (PAC)](https://www.kennedykrieger.org/physiologic-metabolic-anatomic-biomarkers/research/collaborative-projects/preclinical-alzheimer-s-disease-ad-consortium)
- [Open Access Series of Imaging Studies 3 (OASIS-3)](https://sites.wustl.edu/oasisbrains/)

## Scripts

Scripts are organized into the following subdirectories:

- 0_MergeTables: scripts for tidying and merging tabular data from each dataset
- 1_SubjectSelection: scripts for identifying stable and progressor subjects to train SVM models
- 2_Preprocessing: scripts for processing raw PET into regional SUVRs and raw MRI into regional volumes
- 3_TrainTestModel: scripts for training and evaluating SVM models
- 4_FeatureImportance: scripts for computing Haufe-corrected linear SVM feature importance
- 5_A4RetrospectiveAnalysis: scripts for applying trained models on the A4 cohort in the retrospective cohort enrichment experiment
- 6_TablesAndFigures: scripts to generate figures and tables

## Cite

> In revision at *Neurobiology of Aging*
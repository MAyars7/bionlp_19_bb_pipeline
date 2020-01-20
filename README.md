# bionlp_19_bb_pipeline
Pipeline for completing the BioNLP 2019 bacterial biotope NER/norm task

The BioNLP bacterial biotope shared task is to extract bacterial organisms, habitats, and medical phenotypes from Pubmed articles and link them by association.

The Shared Task description is here:
https://drive.google.com/file/d/1G0po_xlRjQCZ-qxuA_4PLdipXU6rtYTp/view

# NER/Norm

The first part of the shared task is Named Entity Recognition (NER) and Normalization.

The goal of this step is to correctly identify and classify words or multi-word phrases in biomedical texts that correspond to the entities of interest in this task: bacterial species, bacterial habitats, medical phenotypes, and geographical locations.

Scripts to generate resources:
---
Entities: \
-generate_bacteria_taxid_dict.py (bacteria) \
-extract_obo_category_nodes.py (habitat, phenotype) \
\
Training data: \
-easy_pubmed_batch_downloads.R (pubmed articles) 

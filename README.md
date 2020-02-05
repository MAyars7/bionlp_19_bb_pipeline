# bionlp_19_bb_pipeline
WIP pipeline for completing the BioNLP 2019 bacterial biotope NER/norm task

The BioNLP bacterial biotope shared task is to extract bacterial organisms, habitats, and medical phenotypes from Pubmed articles and link them by association.

The Shared Task description is here:
https://drive.google.com/file/d/1G0po_xlRjQCZ-qxuA_4PLdipXU6rtYTp/view

# NER/Norm

The first part of the shared task is Named Entity Recognition (NER) and Normalization.

The goal of this step is to correctly identify and classify words or multi-word phrases in biomedical texts that correspond to the entities of interest in this task: bacterial species, bacterial habitats, medical phenotypes, and geographical locations.

The rough steps to accomplish this are:
  1.  Generate lists of the entities of interest from NCBI (bacteria) and the provided .obo resource (habitats, phenotypes).
  2.  Annotate these entities in biomedical texts.
  3.  Train a model on the annotated texts to recognize new entities used in similar contexts to those in the original lists (NER).
  4.  Create general rules to flexibly link clusters of synonymous entities to a single identity (normalization).




Scripts to generate entity lists:
---
-generate_bacteria_taxid_dict.py (bacteria) \
-extract_obo_category_nodes.py (habitat, phenotype) \


Scripts to obtain and annotate biomedical texts: 
---
-easy_pubmed_batch_downloads.R

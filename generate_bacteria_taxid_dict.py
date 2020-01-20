import pandas as pd
import re
import json
import argparse

from ete3 import NCBITaxa

ncbi = NCBITaxa()

parser = argparse.ArgumentParser(
    description="Extract all terms linked to specified concept in obo format file and dump them to a json file.")
parser.add_argument('ncbi_lineage_file',
                    help='Path to ncbi taxdump lineage file.  Each line is a taxonomic lineage.  Ex: ')
parser.add_argument('out_file', help='Path to json to dump microorganism taxid dict')
parser.add_argument('--stopwords_file',
                    help='Stopwords from bioNLP task description.  Species containing any of these words will be filtered out.')
parser.add_argument('--genera_filter_file',
                    help='List of genera to filter to.  If present, output dict will only contain species of these genera.')
args = parser.parse_args()

"""
This script is a component of the BioNLP bacterial biotope named entity recognition/normalization step.

This script generates a dictionary mapping NCBI organism names to their taxonomic identifiers for use in training NLP models to identify, normalize, and link
literature uses of bacterial organisms to their identity.

Shared Task specifications:

https://drive.google.com/file/d/1G0po_xlRjQCZ-qxuA_4PLdipXU6rtYTp/view

To-do:
    -Refine 'corrected name' to account for archaic/redundant/synonymous names.
    -Examine whether names at unnamed NCBI ranks ('No rank 17', 'group', etc.) are used in literature.
    -Extend dict entries to include other lineage ranks.
    -Make additional manual entries from task specifications.
"""

def check_taxon_name_legitimacy(taxon, stopwords, verbose=False):
    """
    Check a taxon name for usability, defined as:
        -Is a string
        -Does not contain a stopword

    To-do:
        Add other conditions for failing taxon names:
            -Rank dependent length (Genus: 1 word, species: 2 words) ?
            -Combinations of 'soft' stop words ?

    :param taxon (string): organism name
    :param stopwords (list): strings that must not be included in taxon
    :param verbose (bool): flag to print reasons for a name failing
    :return: legitimacy (bool): Flag to determine if an organism name is usable
    """

    legitimacy = False

    if type(taxon) == str:
        # if len(taxon.split(' ')) > 1:
        if not any(word in stopwords for word in taxon.lower().split(' ')):
            legitimacy = True
        elif verbose == True:
            print("Stopword in taxon: %s" % taxon)

    return legitimacy

    """
    Iterate through provided names and determine the most specific (e.g. subspecies vs. species), usable name.

    Usability is defined as the check_taxon_name_legitimacy function returning True.

    Input:
        taxid = int, taxonomic identifier
        taxons = list, taxonomic names
    Return:
        (taxid, taxon) = tuple, where taxon is the lastmost usable list entry, or None if no names are usable.
    """

def rightmost_usable_name(taxid, taxons, stopwords):
    """
    From a list of organism names ordered by rank, identify the rightmost usable name.

    :param taxid (int): NCBI taxonomic identifier
    :param taxons (list): NCBI organism names ordered by specificity
    :param stopwords (list): Words that must not be included in a name
    :return (tuple): rightmost usable taxon if there is one, otherwise None
    """
    for taxon in taxons:
        if check_taxon_name_legitimacy(taxon):
            return (taxon, taxid)
    return None

def generate_dict_entry(taxid, varietas, subspecies, species, genus, stopwords):
    """
    This function generates dict entries to be added to the global_microorganism_taxid_dict.

    Provided names for an organism are tested at right-to-left lineage ranks, stopping at the first usable name.

    The idea being that taxids assigned to specific but unusable names will be assigned to higher ranks rather than discarded.

    Ex:
        uncultured Blastocatellia bacterium (2496907) is unusable because it contains the stopword 'uncultured'.  It will be assigned to "Blastocatellia"
        instead, because that is the rightmost usable name in the lineage.

    For species-level names, an additional entry is added to capture the common practice in literature of abbreviating a genus to its first letter

    Ex:
        Staphylococcus aureus -> S. aureus

    :param taxid: int, NCBI taxonomic identifier
    :param varietas: string, varietas rank NCBI organism name
    :param subspecies: string, subspecies rank NCBI organism name
    :param species: string, species rank NCBI organism name
    :param genus: string, genus rank NCBI organism name
    :param stopwords: list of strings, organism names containing any of these strings will be excluded
    """
    global microorganism_taxid_dict

    rank = ncbi.get_rank([taxid]).get(taxid)
    if rank == 'genus':
        if check_taxon_name_legitimacy(genus, stopwords):
            microorganism_taxid_dict[genus] = {
                'taxid': taxid,
                'corrected_name': genus,
                'rank': rank
            }
    elif rank == 'species':
        if check_taxon_name_legitimacy(species, stopwords):

            if species.split(' ')[-1] == 'sp.':
                if check_taxon_name_legitimacy(genus, stopwords):
                    microorganism_taxid_dict[species] = {
                        'taxid': taxid,
                        'corrected_name': genus,
                        'rank': rank
                    }
            else:
                microorganism_taxid_dict[species] = {
                    'taxid': taxid,
                    'corrected_name': species,
                    'rank': rank
                }
                if 'sp.' not in species:
                    abbrev_prefix = '%s.' % species[0]
                    abbrev_species_name = '%s %s' % (abbrev_prefix, ' '.join(species.split(' ')[1:]))
                    microorganism_taxid_dict[abbrev_species_name] = {
                        'taxid': taxid,
                        'corrected_name': species,
                        'rank': rank
                    }
    else:
        for taxon in [varietas, subspecies]:
            if check_taxon_name_legitimacy(taxon, stopwords):
                microorganism_taxid_dict[taxon] = {
                    'taxid': taxid,
                    'corrected_name': taxon,
                    'rank': rank
                }

def generate_truncated_dict_entry(taxid, species, genus, genera_filter_list, stopwords):
    """
    This function limits entries to species-level taxa (only) in a
    select set of genera with pathogenic members.  Genus-abbreviated entries are not included.

    taxid: int, NCBI identifier
    species: str, species name
    genus: str, genus name
    genera_filter_list: list of str, genera
    """
    global microorganism_taxid_dict

    rank = ncbi.get_rank([taxid]).get(taxid)

    if rank == 'species':
        if genus in genera_filter_list:
            if check_taxon_name_legitimacy(species, stopwords):
                microorganism_taxid_dict[species] = {
                    'taxid': taxid,
                    'rank': rank
                }

if __name__ == "__main__":
    """
    Read NCBI lineages from file and filter them to bacterial lineages with taxids belonging to the BioNLP-BB-norm specified list
    of usable taxids.
    """
    print("Generating dataframe of bacterial candidates from NCBI lineages and BioNLP taxids...")
    lineage_df = pd.read_csv(args.ncbi_lineage_file)
    bio_nlp_taxids_path = './resources/BioNLP-OST-2019_BB-norm_Microorganism-ids.txt'

    valid_taxids = [line.rstrip('\n') for line in open(bio_nlp_taxids_path)]
    valid_microorganisms_df = lineage_df[lineage_df['tax_id'].isin(valid_taxids)]
    valid_bacteria_df = valid_microorganisms_df[valid_microorganisms_df['superkingdom'] == 'Bacteria']

    print("Total number of NCBI lineages = %d" % len(lineage_df))
    print("Lineages with valid microorganism taxids = %d" % len(valid_microorganisms_df))
    print("Lineages with superkingdom bacteria and valid microorganism taxids = %d" % len(valid_bacteria_df))

    # If a stopwords file has been specified, read stopwords into list
    stopwords = []
    if args.stopwords_file:
        stopwords = [line.rstrip('\n') for line in open(args.stopwords_file)]
    """
    This will be used as a global dict that entries are added to as they are generated.

    Once completed, it will be written to out_file.
    """
    microorganism_taxid_dict = {}

    """
    If a list of genera to filter to is specified, add only species-level entries to microorganism_taxid_dict for lineages in one of these genera.

    Otherwise, for each lineage, add an entry for the rightmost, usable name (of genus, species, subspecies, varietas).  
    """
    filter_genera = []
    if args.genera_filter_file:

        dict_type = 'truncated'
        filter_genera = [line.rstrip('\n') for line in open(args.genera_filter_file)]
        x = valid_bacteria_df.apply(
            lambda row: generate_truncated_dict_entry(row['tax_id'], row['species'], row['genus'], filter_genera,
                                                      stopwords), axis=1)
    else:
        dict_type = 'full'
        x = valid_bacteria_df.apply(
            lambda row: generate_dict_entry(row['tax_id'], row['varietas'], row['subspecies'], row['species'],
                                            row['genus'], stopwords), axis=1)

    """
    Manual additions:

        These are included in the task specifications.

        “low G+C gram-positive bacteria”, synonym of Firmicutes
        “high G+C gram-positive bacteria”, synonym of Actinobacteria 
    """
    microorganism_taxid_dict['low G+C gram-positive bacteria'] = {

        'taxid': ncbi.get_name_translator(['Firmicutes']).get('Firmicutes'),
        'corrected_name': 'Firmicutes',
        'rank': 'phylum'
    }
    microorganism_taxid_dict['high G+C gram-positive bacteria'] = {
        'taxid': ncbi.get_name_translator(['Actinobacteria']).get('Actinobacteria'),
        'corrected_name': 'Actinobacteria',
        'rank': 'phylum'
    }

    print("Number of entries in %s microorganism taxid dict: %d" % (dict_type, len(microorganism_taxid_dict)))
    print("Writing microorganism taxid dict to %s..." % args.out_file)
    with open(args.out_file, 'w') as f:
        json.dump(microorganism_taxid_dict, f)

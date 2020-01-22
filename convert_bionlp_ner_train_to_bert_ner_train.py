import numpy as np
import pandas as pd
import glob
import scispacy
import spacy
import re
import argparse

parser = argparse.ArgumentParser(description="")
parser.add_argument('bionlp_train_dir', help='Directory containing BioNLP NER train files.')
parser.add_argument('bert_train_outfile', help='Path to write BERT NER train file to.')
args = parser.parse_args()

def extract_and_merge_title_and_abstract(bionlp_lines):
    """
    Given lines from a BioNLP .a1 format file, extract title and abstract (if present) and merge them into a single string.

    :param bionlp_lines (list): lines from a BioNLP .a1 format file
    :return merged_title_abstr (str): merged title (if present) and abstract (if present)
    """
    title_pattern = re.compile('\tTitle ')
    abstr_pattern = re.compile('\tParagraph ')

    title_lines = [i.split('\t')[2].strip() for i in list(filter(title_pattern.search, bionlp_lines))]
    abstr_lines = [i.split('\t')[2].strip() for i in list(filter(abstr_pattern.search, bionlp_lines))]
    merged_title_abstr = ' '.join(title_lines + abstr_lines)

    return merged_title_abstr

def extract_ent_lines(bionlp_lines):
    """
    Extract entity lines which specify the text, label, start position, and end position of entities within the text (1 per line).

    :param bionlp_lines (list): lines from a BioNLP .a1 format file
    :return ent_lines (list): lines specifying a specific entity, rather than a body of text
    """

    #Entity lines are lines that do not contain "Title" or "Paragraph" label
    ent_line_pattern = re.compile(r'^((?!(Title|Paragraph)).)*$')

    ent_lines = [i.strip() for i in list(filter(ent_line_pattern.search, bionlp_lines))]

    return ent_lines

def get_sentence_break_indices(doc_text, nlp):
    """
    Using a scispacy model (trained on biomedical texts), identify the position of breaks between sentences in a text passage.

    :param doc_text (str): biomedical text passage.
    :param nlp (spacy model): model to segment sentences with.
    :return sentence_break_indices (list): positions in the text string of spaces between sentences.
    """

    doc = nlp(doc_text)
    sentence_break_indices = []
    idx = 0
    for sentence in doc.sents:
        idx += len(str(sentence)) + 1
        sentence_break_indices.append(idx)

    return sentence_break_indices

def convert_bionlp_abstract_to_bert_train_format(passage_text, ent_lines, sentence_break_indices):
    """
    Convert the BioNLP NER annotation format into BERT NER annotation format.

    Ex:
        reported	O
        that	O
        β-glucan	B-PHE
        produced	I-PHE
        by	O

    :param passage_text (str): Biomedical text passage
    :param ent_lines (list): Lines specifying an entity annotation (from extract_ent_lines function)
    :param sentence_break_indices (list): Index positions of spaces between sentences in the passage_text string (from get_sentence_break_indices)
    :return out_lines (list): Lines in BERT NER annotation format
    """

    #Conversion of BioNLP labels
    label_dict = {
        'Microorganism': 'MORG',
        'Habitat': 'HAB',
        'Phenotype': 'PHE'
    }

    ent_dict = {}

    """
    Parse ent_lines into ent_dict, 1 entry per token: {
        start index of token in passage_text str: 
            { 
             'end' : end index of token in passage_text string,
             'token' : token (word or symbol),
             'label' : BERT label from label dict, preceded by B- (first token) or I- (second or later token)
            }
        }            
    """
    for ent_line in ent_lines:
        label_start_end = ent_line.split('\t')[1].split(' ')
        label = label_start_end[0]
        start = label_start_end[1]
        end = label_start_end[-1]
        bert_label = label_dict[label]

        start = int(start)
        end = int(end)

        ent_tokens = ent_line.split('\t')[2].split(' ')

        #first token of entity
        ent_dict[start] = {
            'end': end,
            'token': ent_tokens[0],
            'label': 'B-%s' % bert_label
        }
        idx = start + len(ent_tokens[0]) + 1

        #Tokens after first token of entity
        for token in ent_tokens[1:]:
            ent_dict[idx] = {
                'end': idx + len(token),
                'token': token,
                'label': 'I-%s' % bert_label
            }
            idx += len(token) + 1

    out_lines = []
    abstract_idx = 0

    """
    Iterate through words of passage_text string, generating 1 line per word.
    
    Non-entity tokens are labeled 'O', entities are labeled 'B-<label>' (first token of entity) or 'I-<label>' (second or later token of entity)
    """
    for word in passage_text.split(' '):
        if word:
            ent = ent_dict.get(abstract_idx)
            end = abstract_idx + len(word) + 1

            if ent:
                out_label = ent['label']
            else:
                out_label = 'O'

            """
            If a word ends in a period and is the last word of a sentence, generate two lines: 1 for the word, 1 for the period.
            
            Otherwise, generate 1 line for the word.
            """
            if word[-1] == '.' and end in sentence_break_indices:
                out_lines.append('%s\t%s\n' % (word[:-1], out_label))
                out_lines.append('%s\tO\n' % (word[-1]))
            else:
                out_lines.append('%s\t%s\n' % (word, out_label))

            abstract_idx += len(word) + 1
            if abstract_idx in sentence_break_indices:
                out_lines.append('\n')

    return out_lines

def generate_article_train_lines(article_text, bionlp_lines, nlp):
    """
    Convert a single BioNLP annotated title+abstract or body passage into BERT NER training format.

    :param article_text (str): Biomedical text passage
    :param bionlp_lines (list): BioNLP annotation lines
    :param nlp (spacy model): model to segment sentences with
    :return article_train_lines (list): BERT NER annotation lines
    """

    ent_lines = extract_ent_lines(bionlp_lines)
    sentence_break_indices = get_sentence_break_indices(article_text, nlp)
    article_train_lines = convert_bionlp_abstract_to_bert_train_format(article_text, ent_lines, sentence_break_indices)

    return article_train_lines

if __name__ == "__main__":
    nlp = spacy.load('en_core_sci_md')

    out_file = args.bert_train_outfile

    """
    BioNLP NER annotation files exist in 2 formats:
        
        Format 1:
            <BB-norm-PMID> - Single .a1 file contains labeled Title and/or Abstract followed by lines with entity annotations.
        Format 2:
            <BB-norm-F-PMID-NNN> - Pairs of txt files (containing the text passage) and .a1 files (containing entity annotations)
    
    Regex patterns are used to differentiate the two formats, and pair up the text and a1 files for Format 2.
    """
    bionlp_train_files = glob.glob('%s*a1' % args.bionlp_train_dir)
    text_files = glob.glob('%s*txt' % args.bionlp_train_dir)

    pm_abstr_pattern = re.compile(r'BB-norm\-\w+.a1')
    pm_passage_a1_pattern = re.compile(r'BB-norm-F\-\w+\-\w{3}.a1')
    pm_passage_txt_pattern = re.compile(r'BB-norm-F\-\w+\-\w{3}.txt')

    abstr_files = list(filter(pm_abstr_pattern.search, bionlp_train_files))
    passage_a1_files = sorted(list(filter(pm_passage_a1_pattern.search, bionlp_train_files)))
    passage_txt_files = sorted(list(filter(pm_passage_txt_pattern.search, text_files)))

    bert_train_lines = []
    bert_train_lines_2 = []

    usable_files = 0
    unusable_files = []

    # Title/abstract files, Format 1
    for bionlp_train_file in abstr_files:

        with open(bionlp_train_file) as f:
            bionlp_lines = f.readlines()

        if len(bionlp_lines) > 1:
            usable_files += 1

            merged_title_abstract = extract_and_merge_title_and_abstract(bionlp_lines)
            bert_train_lines += generate_article_train_lines(merged_title_abstract, bionlp_lines, nlp)

        else:
            unusable_files.append(bionlp_train_file)

    # Passage files, Format 2
    passage_file_tuples = list(zip(sorted(passage_a1_files), sorted(passage_txt_files)))

    for passage_file_tuple in passage_file_tuples:

        a1_file, txt_file = passage_file_tuple

        with open(a1_file) as f:
            bionlp_lines = f.readlines()

        with open(txt_file) as f:
            passage_text = ' '.join([i.strip() for i in f.readlines()])

        bert_train_lines += generate_article_train_lines(passage_text, bionlp_lines, nlp)

    print("Directory contains %d Title/Abstract files." % len(abstr_files))
    print("Directory contains %d body passage files." % len(passage_file_tuples))

    with open(out_file, 'w') as f:
        for line in bert_train_lines:
            f.write(line)
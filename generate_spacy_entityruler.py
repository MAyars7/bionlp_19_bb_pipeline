import json
import spacy
import en_core_web_sm
from spacy.pipeline import EntityRuler
from spacy import displacy

import itertools
import re

import networkx

if __name__ == '__main__':

    #Load Spacy english model and create an EntityRuler object
    nlp = en_core_web_sm.load()
    ruler = EntityRuler(nlp)

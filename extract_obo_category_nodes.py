import json
import itertools
import re
import networkx
import argparse

"""
This script is a component of the BioNLP bacterial biotope named entity recognition/normalization step.

Shared Task specifications:
https://drive.google.com/file/d/1G0po_xlRjQCZ-qxuA_4PLdipXU6rtYTp/view

This script contains functions to extract nodes linked to a specified concept/category from an obo file into a graph object,
convert the graph into a dict, and dump the dict to a json file. 

To-do:
    -Add ability to extract by a list of categories
    -Add functions to visualize paths between graph nodes
    -Add handling for non-exact synonyms
"""
parser = argparse.ArgumentParser(description="Extract all terms linked to specified concept in obo format file and dump them to a json file.")
parser.add_argument('obo_file_path', help='Path to .obo format file')
parser.add_argument('category', help='Name of category/parent/concept node.  All extracted nodes will be linked to this node.  Ex: "microbial habitats"')
parser.add_argument('out_file', help='Name of json file to dump dict of extracted nodes to.')
args = parser.parse_args()

def parse_obo_entry(entry_lines):
    """
    Convert entry lines from obo file into a dict entry.
    :param entry_lines (list): strings extracted from .obo file for a
    :return: entry (dict), each regex match in the term is a key/value entry
    """
    entry_type = next(entry_lines)
    entry_lines = list(entry_lines)
    entry = dict()

    # Regex pattern for parsing obo format from https://github.com/cmungall/obo/blob/master/obo/read.py
    obo_entry_line_pattern = re.compile(
        r'^(?P<key>.+?): *(?P<value>.+?) ?(?P<trailing_modifier>(?<!\\)\{.*?(?<!\\)\})? ?(?P<comment>(?<!\\)!.*?)?$')

    #Verify entry is a Term, regex match key and value, and add it to dict
    if entry_type.startswith('[Term]'):
        for line in entry_lines:
            if line.startswith('!'):
                continue
            regex_match = re.match(obo_entry_line_pattern, line)
            key = regex_match.group('key')
            value = regex_match.group('value')
            entry.setdefault(key, []).append(value)

    return entry

def parse_obo_file_to_graph(obo_path):
    """
    Parse obo file into a networkx MultiDiGraph.

    Regex pattern for parsing obo terms and some parts of this function are copied from:
        https://github.com/cmungall/obo/blob/master/obo/read.py

    :param obo_path (string): path to .obo file.
    :return: graph (MultiDiGraph): graph containing linked entries from .obo file
    """
    graph = networkx.MultiDiGraph()
    edge_tuples = []

    with open(obo_path) as f:
        obo_lines = f.read().splitlines()

    entries = itertools.groupby(obo_lines, lambda line: line.strip() == '')
    parsed_terms = []
    for is_blank, entry_lines in entries:
        if is_blank:
            continue
        parsed_term = parse_obo_entry(entry_lines)
        parsed_terms.append(parsed_term)

    for term in parsed_terms:
        if term:
            term_id = term.pop('id')[0]
            graph.add_node(term_id, attr_dict=term)
            for target_term in term.pop('is_a', []):
                edge_tuple = term_id, 'is_a', target_term
                edge_tuples.append(edge_tuple)

    for origin_term, typedef, dest_term in edge_tuples:
        graph.add_edge(origin_term, dest_term, key=typedef)

    return graph

def get_nodes_by_category(graph, category, convert_to_names=False):
    """
    Returns a list of nodes linked to a specific category from a MultiDiGraph object.

    :param graph (MultiDiGraph): MultiDiGraph from parse_obo_file_to_graph function
    :param category (string): name of category/parent node to get linked nodes from, ex: 'microbial habitat'
    :param convert_to_names (bool): If True, returns a list of name strings instead of node ids
    :return: category nodes (list): node ids or node names for specific category
    """
    category_nodes = []

    for source_node in graph.nodes():
        for i in networkx.algorithms.traversal.edgedfs.edge_dfs(graph, source=source_node):
            node_name = (graph.nodes.get(i[1]).get('attr_dict').get('name')[0])
            if node_name == category:
                category_nodes.append(source_node)

    # Convert category nodes from list of identifiers (OBT:002837) to string names
    if convert_to_names:
        category_nodes = [graph.nodes.get(i).get('attr_dict').get('name')[0] for i in category_nodes]

    return category_nodes

def convert_node_list_to_dict(graph, node_list, category, path_to_write=''):
    """
    :param graph (MultiDiGraph): from parse_obo_file_to_graph function
    :param node_list (list): string node ids belonging to graph
    :param category (string): name of node all nodes in list are linked to
    :param path_to_write (string): path to write dict to
    :return: node_dict (dict): {node_name : { 'id' : node id, 'entity_class' : category argument}}
    """
    node_dict = {}

    obo_entry_synonym_pattern = re.compile(
        r'(\"[ \w|\W+]+\") ([A-Z]+) (\[[a-zA-Z]*\:*[0-9]*\])'
    )

    # Parse nodes into dict with ID, name, and synonyms.  For now, only use synonyms of type EXACT.
    for node_id in node_list:
        node = graph.nodes.get(node_id)

        # Obtain node attribute fields.  Node ID (key) and node name are required, synonyms is optional.
        assert node, "Node %s does not exist." % node_id

        node_attr_dict = node.get('attr_dict')
        assert node_attr_dict, "Node %s has no attached attribute dict." % node

        node_names = node_attr_dict.get('name')
        assert node_names, "Node %s has no name field." % node

        # Possible for nodes to have more than 1 name?
        node_name = node_names[0]

        node_dict[node_name] = {
            'id': node_id,
            'entity_class': category
        }
        node_synonyms = []
        attr_dict_node_synonyms = node_attr_dict.get('synonym')
        if attr_dict_node_synonyms:
            for attr_synonym in attr_dict_node_synonyms:
                attr_synonym = attr_synonym.replace('-', ' ')
                re_match_syn = re.match(obo_entry_synonym_pattern, attr_synonym)
                re_match_syn_groups = re_match_syn.groups()

                attr_synonym_name = re_match_syn_groups[0].replace('"', '')
                attr_synonym_type = re_match_syn_groups[1]
                attr_synonym_id = re_match_syn_groups[2]

                if attr_synonym_type == 'EXACT':
                    node_dict[attr_synonym_name] = {
                        'id': node_id,
                        'entity_class': category
                    }
    return node_dict

if __name__ == "__main__":

    print("Generating MultiDiGraph from obo file...")
    graph = parse_obo_file_to_graph(args.obo_file_path)

    category = args.category.replace('_', ' ')
    print("Extracting nodes linked to category %s..." % category)
    category_nodes = get_nodes_by_category(graph, category)

    print("Converting extracted nodes into dict...")
    category_nodes_dict = convert_node_list_to_dict(graph, category_nodes, args.category)

    print("Writing category nodes dict to %s..." % args.out_file)
    with open(args.out_file, 'w') as f:
        json.dump(category_nodes_dict, f)

import pandas as pd
import numpy as np
import json
from scipy import stats


level = [
    'superkingdom', 'phylum', 'class',
    'order', 'family', 'genus', 'species'
]

def match_peptide_taxonomy(quantity, taxonomy, equate_il = True):
    taxonomy = taxonomy[
        [x + '_name' for x in level] + ['peptide']
    ]
    taxonomy = taxonomy.rename(columns={x + '_name': x for x in level})

    quantity = quantity[
        ['PEP.StrippedSequence'] + quantity.columns[
            quantity.columns.str.contains('.PEP.Quantity', regex=False)
        ].tolist()
    ]
    quantity = quantity.rename(
        columns={
            x: x.replace('.PEP.Quantity', '.Quantity')
            for x in quantity.columns[1:]
        }
    )
    quantity.rename(columns={'PEP.StrippedSequence': 'peptide'}, inplace=True)

    if equate_il:
        taxonomy.loc[:, 'peptide'] = taxonomy['peptide'].str.replace('I', 'L')
        quantity.loc[:, 'peptide'] = quantity['peptide'].str.replace('I', 'L')

    data = pd.merge(
        taxonomy, quantity,
        on='peptide', how='inner'
    )
    return data


def clean_taxonomy(data, min_peptide_count = 2):
    for i in range(len(level) - 1, 0, -1):
        count = data[level[i]].value_counts()
        taxa = count[count < min_peptide_count].index
        data.loc[data[level[i]].isin(taxa), level[i]] = np.NaN

    for i in range(len(level) - 1, 0, -1):
        data.loc[~data[level[i]].isnull(), level[i - 1]] = \
            data.loc[~data[level[i]].isnull(), level[i - 1]] \
                .fillna('-')

    data.dropna(subset=[level[0]], axis=0, inplace=True)
    data.fillna({l: '' for l in level}, inplace=True)
    return data


def build_tree(data):
    data = data.groupby(by=level).sum()

    tree = {'name': 'root', 'children': []}

    for key, value in data.iterrows():
        depth_cursor = tree

        for depth, lvl in enumerate(level):
            index = None

            for i, child in enumerate(depth_cursor['children']):
                if key[depth] == child['name']:
                    index = i

            if index is None:
                depth_cursor['children'] \
                    .append({'name': key[depth], 'children': []})
                index = len(depth_cursor['children']) - 1

            depth_cursor = depth_cursor['children'][index]

            if depth == len(level) - 1:
                depth_cursor['value'] = value

    return tree


def calculate_sum(tree):
    value = tree.get('value', None)
    if value is not None:
        return value

    children = tree.get('children', None)
    if not children:
        return None

    for child in children:
        child_value = calculate_sum(child)
        if child_value is not None:
            if value is None:
                value = child_value
            else:
                value = value + child_value

    tree['value'] = value
    return value


def cut_tree(tree, depth):
    children = tree.get('children', None)

    if depth > 0:
        if children:
            for child in children:
                cut_tree(child, depth - 1)
    else:
        if children:
            children = []
            tree['children'] = children

    if tree['name'] == '-':
        if not children or \
            not any(map(lambda c: c['name'] != '', children)):
            tree['name'] = ''


def calculate_fc_pvalue(tree, group, method = 'ttest', **test_args):
    index1 = group[group == 1].index
    index2 = group[group == 2].index

    value1 = tree['value'][index1]
    value2 = tree['value'][index2]

    mean1 = value1.mean()
    mean2 = value2.mean()
    fc = np.divide(mean2, mean1)

    
    if method == 'ttest':
        test = stats.ttest_ind(value1, value2, **test_args)
    elif method == 'mannwhitneyu':
        test = stats.mannwhitneyu(
            values[j], values[i],
            **test_args
        )
    else:
        raise ValueError('invalid method')
    pvalue = test.pvalue

    tree['data'] = {
        'value1': value1.tolist(),
        'value2': value2.tolist(),
        'mean1': float(mean1),
        'mean2': float(mean2),
        'fc': float(fc),
        'pvalue': float(pvalue)
    }

    children = tree.get('children', None)
    if children:
        for child in children:
            calculate_fc_pvalue(child, group, **test_args)

    return tree['data']


def finalize_tree_for_output(tree):
    children = tree.get('children', None)

    if children:
        for child in children:
            finalize_tree_for_output(child)

    del tree['value']


def save_tree_json(tree, file, jsonp=False):
    with open(file, 'w') as f:
        if jsonp:
            if not isinstance(jsonp, str):
                jsonp = 'getData'
            f.write(jsonp + '(')

        json.dump(tree, f, skipkeys=True)

        if jsonp:
            f.write(')')


def convert_tree_to_dataframe(tree):
    def _to_dataframe(node, depth):
        if depth == 0:
            lvl = 'root'
        else:
            lvl = level[depth - 1]

        data = node['data']

        result = pd.DataFrame.from_records([
            pd.Series([lvl, node['name']], index=['level', 'name']).append([
                pd.Series(
                    data[name],
                    index=[name + '.' + str(i + 1) for i, x in enumerate(data[name])]
                        if isinstance(data[name], list)
                        else [name]
                )
                for name in ['value1', 'value2', 'mean1', 'mean2', 'fc', 'pvalue']
            ])
        ])

        children = node.get('children', None)
        if children:
            children_result = [
                _to_dataframe(child, depth=depth + 1)
                for child in children
                if child['name'] != '' and child['name'] is not None
            ]
            if len(children_result) > 0:
                result = result.append(children_result, ignore_index=True)

        return result

    return _to_dataframe(tree, 0)



if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Generate taxonomy tree for visualization'
    )
    parser.add_argument(
        '--quantity', required=True,
        help='input peptide report file'
    )
    parser.add_argument(
        '--taxonomy', required=True,
        help='input unipept taxonomy file'
    )
    parser.add_argument(
        '--out', default='taxonomy_tree.json',
        help='output taxonomy tree file (default: %(default)s)'
    )
    parser.add_argument(
        '--group', required=True,
        help='group for each sample in the peptide report file, 1 or 2 separated with ","'
    )
    parser.add_argument(
        '--min_peptide_count', type=int, default=2,
        help='the minimum number of peptides required for a taxon (default: %(default)s)'
    )
    parser.add_argument(
        '--depth', type=int, default=7,
        help='the depth of taxonomy tree (default: %(default)s)'
    )
    parser.add_argument(
        '--jsonp', default='getData',
        help='JSONP function name in the output file (default: %(default)s)'
    )
    parser.add_argument(
        '--out_csv', default='taxonomy_quantity.csv',
        help='output taxonomy quantity table file (default: %(default)s)'
    )


    args = parser.parse_args()
    quantity_file = args.quantity
    taxonomy_file = args.taxonomy
    out_file = args.out
    min_peptide_count = args.min_peptide_count
    depth = args.depth
    group = args.group
    jsonp = args.jsonp
    out_csv_file = args.out_csv

    quantity = pd.read_csv(quantity_file)
    taxonomy = pd.read_csv(taxonomy_file)

    data = match_peptide_taxonomy(quantity, taxonomy)
    data = clean_taxonomy(data, min_peptide_count = min_peptide_count)

    tree = build_tree(data)
    calculate_sum(tree)
    cut_tree(tree, depth)

    group = pd.Series(
        map(int, group.split(',')),
        index=data.columns[len(level) + 1:]
    )
    calculate_fc_pvalue(tree, group)

    finalize_tree_for_output(tree)
    save_tree_json(tree, out_file, jsonp=jsonp)

    if out_csv_file is not None:
        taxonomy_quantity = convert_tree_to_dataframe(tree)
        taxonomy_quantity.to_csv(out_csv_file, index=False)



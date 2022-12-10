import pandas as pd
import numpy as np
import itertools
from scipy import stats
from statsmodels.stats import multitest

from taxonomy_tree import *

level = [
    'superkingdom', 'phylum', 'class',
    'order', 'family', 'genus', 'species'
]


def enumerate_tree(tree):
    stack = []
    stack.append(tree)
    while stack:
        node = stack.pop()
        yield node            
        for child in node['children']:
            stack.append(child)


def calculate_percentage(tree):
    root_value = tree['value']
    for node in enumerate_tree(tree):
        value = node.get('value', None)
        if value is not None:
            node['value'] = value / root_value * 100
    return tree


def calculate_fc_pvalue_multigroups(tree, group):
    values = [
        tree['value'][group[group == g].index]
        for g in sorted(set(group))
    ]
    
    means = [
        value.mean()
        for value in values
    ]
   
    comparisons = list(itertools.combinations(
        range(len(values)), 2
    ))
    
    fold_changes = [
        np.divide(means[j], means[i])
        for i, j in comparisons
    ]

    if np.max(list(map(np.max, values))) == np.min(list(map(np.min, values))):
        kw_pvalue = np.nan
    else:
        kw_pvalue = stats.kruskal(*values).pvalue
    
    mwu_pvalues = [
        stats.mannwhitneyu(
            values[j], values[i],
            alternative='two-sided'
        ).pvalue
        for i, j in comparisons
    ]

    tree['data'] = {
        'value' + str(i + 1): value.tolist()
        for i, value in enumerate(values)
    }    
    tree['data'].update({
        'mean' + str(i + 1): float(mean)
        for i, mean in enumerate(means)
    })
    tree['data'].update({
        'fc' + str(cmp[1] + 1) + 'vs' + str(cmp[0] + 1): float(fc)
        for cmp, fc in zip(comparisons, fold_changes)
    })
    tree['data'].update({
        'pvalue' + str(cmp[1] + 1) + 'vs' + str(cmp[0] + 1): float(fc)
        for cmp, fc in zip(comparisons, mwu_pvalues)
    })

    tree['data'].update({
        'pvalueKW': float(kw_pvalue)
    })
    
    children = tree.get('children', None)
    if children:
        for child in children:
            calculate_fc_pvalue_multigroups(child, group)

    return tree['data']


def adjust_pvalue(tree, method='fdr_bh'):
    for key in filter(
        lambda x: 'pvalue' in x and 'adjustedPValue' not in x, 
        list(tree['data'].keys())
    ):
        raw_pvalues = [
            node['data'][key]
            for node in enumerate_tree(tree)
        ]
        adjusted_pvalues = multitest.multipletests(
            raw_pvalues,
            method=method
        )[1]

        for i, node in enumerate(enumerate_tree(tree)):
            node['data'][key.replace('pvalue', 'adjustedPValue')] = \
                float(adjusted_pvalues[i])

    return tree['data']


def convert_tree_to_dataframe_multigroups(tree):
    keys = tree['data'].keys()
    statnames = [
        filter(lambda k: k.startswith('fc'), keys),
        filter(lambda k: k.startswith('pvalue') and \
            not k.endswith('KW'), keys)
    ]
    if any(filter(lambda x: 'adjustedPValue' in x, keys)):
        statnames.append(filter(lambda k: \
            k.startswith('adjustedPValue') and \
            not k.endswith('KW'), keys))

    value_names = list(itertools.chain(
        filter(lambda k: k.startswith('value'), keys),
        filter(lambda k: k.startswith('mean'), keys),
        (item for pair in zip(*statnames) for item in pair),
        filter(lambda x: x in keys, ['pvalueKW', 'adjustedPValueKW'])
    ))
    
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
                for name in value_names                
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
        '--quantity_type', choices=['pivot', 'long'], default='pivot',
        help='type of input peptide report file (default: %(default)s)'
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
        help='group for each sample in the peptide report file: 1, 2, 3, ... separated with ","'
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
        '--quantify_method', choices=['quantity', 'count'], default='quantity',
        help='method used for quantification (default: %(default)s)'
    )
    parser.add_argument(
        '--abundance', choices=['raw', 'percentage'], default='raw',
        help='whether raw or percentage abundance used for quantification (default: %(default)s)'
    )
    parser.add_argument(
        '--pvalue_adjust_method', choices=['none', 'bonferroni', 'fdr_bh'], default='fdr_bh',
        help='method used for testing and adjustment of pvalues (default: %(default)s)'
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
    pvalue_adjust_method = args.pvalue_adjust_method
    jsonp = args.jsonp
    out_csv_file = args.out_csv
    quantity_type = args.quantity_type
    quantify_method = args.quantify_method
    abundance_type = args.abundance
    
    quantity = pd.read_csv(quantity_file, low_memory=False)
    taxonomy = pd.read_csv(taxonomy_file, low_memory=False)

    if quantity_type == 'long':
        quantity = quantity.pivot(
            index='PEP.StrippedSequence',
            columns='R.FileName',
            values='PEP.Quantity'
        )
        quantity.columns = [
            x + '.PEP.Quantity'
            for x in quantity.columns
        ]
        quantity.reset_index(inplace=True)
        
        quantity.to_csv(quantity_file.replace('.csv', '') + '.pivot.csv', index=False)

    quantity.fillna(
        0,
        # quantity.min(numeric_only=True).min() / 2, 
        inplace=True
    )
    
    if quantify_method == 'count':  
        quantity_columns = quantity.columns.str.contains('.PEP.Quantity', regex=False)    
        quantity.loc[:, quantity_columns] = quantity.loc[:, quantity_columns].applymap(lambda x: int(x > 0))
    
    data = match_peptide_taxonomy(quantity, taxonomy)
    data = clean_taxonomy(data, min_peptide_count = min_peptide_count)
    
    tree = build_tree(data)
    calculate_sum(tree)
    cut_tree(tree, depth)
    
    if abundance_type == 'percentage':
        calculate_percentage(tree)
    
    group = pd.Series(
        map(int, group.split(',')),
        index=data.columns[len(level) + 1:]
    )
    calculate_fc_pvalue_multigroups(tree, group)
    
    if pvalue_adjust_method and pvalue_adjust_method != 'none':
        adjust_pvalue(tree, method=pvalue_adjust_method)

    finalize_tree_for_output(tree)
    save_tree_json(tree, out_file, jsonp=jsonp)

    if out_csv_file is not None:
        taxonomy_quantity = convert_tree_to_dataframe_multigroups(tree)
        taxonomy_quantity.to_csv(out_csv_file, index=False)



import requests
import time
import pandas as pd
import logging


def unipept_pept2lca(peptide, equate_il=False, extra=False, names=False):
    if peptide is None:
        raise ValueError('peptide is None')
    
    url = 'http://api.unipept.ugent.be/api/v1/pept2lca'
    
    data = [
        ('equate_il', 'true' if equate_il else 'false'),
        ('extra', 'true' if extra else 'false'),
        ('names', 'true' if names else 'false')
    ]    
    if isinstance(peptide, str):
        data.append(('input[]', peptide))
    else:
        data += [('input[]', pep) for pep in peptide]
    
    res = requests.post(
        url, 
        headers = {'Accept': 'application/json'},
        data = data
    )
    
    res.raise_for_status()
    return res.json()
    
    
def find_taxonomy(peptide, equate_il=False,
                  batch_size=100, sleep=2,
                  return_generator=False,
                  logger=None):
    def _batch(batch_id, peptide, wait=None):
        if wait is not None and wait > 0:
            time.sleep(wait)
            
        try:
            data = unipept_pept2lca(
                peptide, 
                equate_il=equate_il, 
                extra=True, names=True
            )
            if logger is not None:
                logger.info(
                    '({batch_id}) Taxonomy found for {m} of {n} peptide(s)' \
                    .format(
                        batch_id=batch_id, 
                        m=len(data), n=len(peptide)
                    )
                )
                
        except Exception as e:
            data = None
            if logger is not None:
                logger.error(
                    '({batch_id}) Failed for peptide(s) {peptide}: {err}' \
                    .format(
                        batch_id=batch_id, 
                        peptide=peptide, 
                        err=e.args
                    )
                )
                
        return data
         
    if isinstance(peptide, str):
        peptide = [peptide]
        
    if batch_size is None or batch_size <= 0:
        batch_size = len(peptide)

    result = (
        _batch(
            batch_id='{i}/{n}'.format(
                i=min(i + batch_size, len(peptide)), 
                n=len(peptide)
            ),
            peptide=peptide[i:(i + batch_size)], 
            wait=sleep if i > 0 else None
        )
        for i in range(0, len(peptide), batch_size)
    )
        
    if not return_generator:
        result = list(result)

    return result
    
   

if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Find taxonomy for peptides using Unipept API'
    )
    parser.add_argument(
        '--in', required=True,
        help='input peptide report file'
    )
    parser.add_argument(
        '--out', default='unipept.csv',
        help='output file (default: %(default)s)'
    )
    parser.add_argument(
        '--exclude',
        help='exclude peptide list file'
    )
    parser.add_argument(
        '--skip', type=int, default=0,
        help='the number of peptides skipped from the start of the input file (default: %(default)s)'
    )
    parser.add_argument(
        '--batch_size', type=int, default=100,
        help='the maximum number of peptides queried for each HTTP request (default: %(default)s)'
    )
    parser.add_argument(
        '--sleep', type=float, default=2,
        help='the seconds to wait between two HTTP requests (default: %(default)s)'
    )
    equate_il_group = parser.add_mutually_exclusive_group(required=False)
    equate_il_group.add_argument(
        '--equate_il', 
        dest='equate_il', action='store_true',
        help='isoleucine (I) and leucine (L) are equated when matching tryptic peptides to UniProt entries (default: %(default)s)'
    )
    equate_il_group.add_argument(
        '--no-equate_il', 
        dest='equate_il', action='store_false',
        help='isoleucine (I) and leucine (L) are not equated'
    )
    parser.set_defaults(equate_il=True)
    
    args = parser.parse_args()
    in_file = getattr(args, 'in')
    out_file = args.out
    exclude_file = args.exclude
    equate_il = args.equate_il 
    batch_size = args.batch_size
    sleep = args.sleep
    skip = args.skip
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    peptide = pd.read_csv(in_file)['PEP.StrippedSequence'] \
        .drop_duplicates()
    if skip is not None and skip > 0:
        peptide = peptide[skip:]

    
    if exclude_file is not None:        
        exclude = pd.read_csv(exclude_file)['peptide']
        peptide = peptide[~peptide.isin(exclude)]

    peptide = peptide.tolist()
    
    first = True
    for data in find_taxonomy(
        peptide, equate_il=equate_il,
        batch_size=batch_size, sleep=sleep,
        return_generator=True,
        logger=logging
    ):
        if data is None or len(data) == 0:
            continue
        
        df = pd.DataFrame.from_records(data)
        
        columns = [
            'peptide', 'taxon_id', 'taxon_name', 'taxon_rank',
            'superkingdom_id', 'superkingdom_name',
            'kingdom_id', 'kingdom_name',
            'subkingdom_id', 'subkingdom_name',
            'superphylum_id', 'superphylum_name',
            'phylum_id', 'phylum_name',
            'subphylum_id', 'subphylum_name',
            'superclass_id', 'superclass_name',
            'class_id', 'class_name',
            'subclass_id', 'subclass_name',
            'infraclass_id', 'infraclass_name',
            'superorder_id', 'superorder_name',
            'order_id', 'order_name',
            'suborder_id', 'suborder_name',
            'infraorder_id', 'infraorder_name',
            'parvorder_id', 'parvorder_name',
            'superfamily_id', 'superfamily_name',
            'family_id', 'family_name',
            'subfamily_id', 'subfamily_name',
            'tribe_id', 'tribe_name',
            'subtribe_id', 'subtribe_name',
            'genus_id', 'genus_name',
            'subgenus_id', 'subgenus_name',
            'species_group_id', 'species_group_name',
            'species_subgroup_id', 'species_subgroup_name',
            'species_id', 'species_name',
            'subspecies_id', 'subspecies_name',
            'varietas_id', 'varietas_name',
            'forma_id', 'forma_name'
        ]
        df = df[columns]

        if first:
            df.to_csv(out_file, index=False)
            first = False
        else:
            df.to_csv(out_file, index=False, mode='a', header=False)
            
            
"""Resolve missing external-cohort symbols without merging unrelated genes."""
from pathlib import Path
import numpy as np
import pandas as pd
from input_utils import canonical_symbols


def external_symbols(var, return_audit=False):
    h = pd.read_csv(Path(__file__).parent / 'reference/hgnc_complete_set.txt', sep='\t', dtype=str).fillna('')
    mapping = h[h.ensembl_gene_id != ''].groupby('ensembl_gene_id').symbol.agg(set).to_dict()
    old = var['gene_symbol'].astype(str).tolist() if 'gene_symbol' in var else var.index.astype(str).tolist()
    result, audit = [], []
    for identifier, symbol in zip(var.index.astype(str), old):
        symbol = symbol.strip()
        key = identifier.split('.')[0]
        if symbol.lower() not in ['', 'nan', 'none', 'na']:
            new, method = symbol, 'existing_symbol'
        elif len(mapping.get(key, set())) == 1:
            new, method = next(iter(mapping[key])), 'unique_HGNC_Ensembl_recovery'
        else:
            new, method = key, 'retain_unresolved_Ensembl_identifier'
        result.append(new)
        audit.append(dict(feature_id=identifier, previous_symbol=symbol, resolved_identifier=new, method=method))
    result = canonical_symbols(result)
    assert all(str(g).strip() and str(g).lower() not in ['nan', 'none', 'na'] for g in result)
    if return_audit:
        frame = pd.DataFrame(audit)
        frame['resolved_identifier'] = result
        return result, frame
    return result

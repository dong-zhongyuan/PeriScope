"""Restore documented EBI transformations before attaching Oxford frequencies.

GRCh38 outcome coordinates remain GRCh38. ST9's paired coordinates are used
only for the exact-variant frequency lookup, never as a generic liftover.
"""
from collections import defaultdict
import numpy as np
import pandas as pd

from .expanded import _chrom

SWAPPED = {2,4,6,8,11,13}
REVERSE = {3,4,7,8,12,13}
SUPPORTED = {1,2,3,4,5,6,7,8,10,11,12,13}
COMP = str.maketrans('ACGT','TGCA')


def normalize_native(rows):
    """Rename the verified Oxford GRCh37 extraction; do not flip any allele."""
    rename = {
        'chromosome': 'chrom_y', 'base_pair_location': 'pos_y',
        'effect_allele': 'ea_y', 'other_allele': 'oa_y', 'beta': 'beta_y',
        'standard_error': 'se_y', 'p_value': 'p_y',
        'neg_log10_p_value': 'neg_log10_p_y', 'variant_id': 'variant_y',
    }
    fields = list(rename) + ['phenotype', 'source_y', 'build_y', 'n_y', 'n_y_status']
    missing = set(fields) - set(rows.columns)
    if missing:
        raise ValueError(f'Missing native outcome fields: {sorted(missing)}')
    if set(rename.values()) & set(rows.columns):
        raise ValueError('Ambiguous mixture of raw-native and normalized columns')
    if not rows.build_y.eq('GRCh37').all():
        raise ValueError('Native Oxford extraction must retain verified GRCh37 coordinates')
    result = rows.loc[:, fields].rename(columns=rename).copy()
    result['chrom_y'] = result.chrom_y.map(_chrom)
    return result


def restore_original(ea,oa,beta,code):
    """Inverse of the documented hm_code strand and allele transformations."""
    if code not in SUPPORTED:
        return ea,oa,np.nan,'unresolved_harmonizer_code'
    if code in SWAPPED: ea,oa,beta=oa,ea,-beta
    if code in REVERSE:
        ea,oa=ea.translate(COMP),oa.translate(COMP)
    return ea,oa,beta,'documented_inverse'


def normalize_ebi(rows, instruments, frequencies):
    positions=defaultdict(set)
    for r in instruments.to_dict('records'):
        positions[(_chrom(r['chrom']),int(r['pos38']),str(r['variant']))].add(int(r['pos37']))
    freq=defaultdict(list)
    for r in frequencies.to_dict('records'):
        freq[(_chrom(r['chr']),int(r['pos37']))].append(r)
    records=[]
    for raw in rows.to_dict('records'):
        code=int(raw['hm_code']) if pd.notna(raw['hm_code']) else -1
        ea,oa,beta,status=restore_original(str(raw['effect_allele']),str(raw['other_allele']),float(raw['beta']),code)
        chrom,pos=_chrom(raw['chromosome']),int(raw['base_pair_location'])
        variants={str(raw[k]) for k in ['variant_id','rsid'] if pd.notna(raw.get(k))}
        mapped=set().union(*(positions.get((chrom,pos,v),set()) for v in variants))
        row=dict(phenotype=raw['phenotype'],build_y='GRCh38',chrom_y=chrom,pos_y=pos,
                 variant_y=raw['variant_id'],ea_y=ea,oa_y=oa,beta_y=beta,se_y=raw['standard_error'],
                 neg_log10_p_y=raw['neg_log10_p_value'],p_y=raw['p_value'],n_y=raw['n_y'],
                 n_y_status=raw['n_y_status'],source_y=raw['source_y'],
                 eaf_y=np.nan,eaf_source=None,eaf_attachment_status='no_exact_paired_coordinate_variant',
                 eaf_reference_cohort='BIG40 release2 genotype frequency resource',
                 hm_code=code,orientation_recovery=status,
                 harmonized_ea=raw['effect_allele'],harmonized_oa=raw['other_allele'],harmonized_beta=raw['beta'])
        if len(mapped)==1 and status=='documented_inverse':
            matches=[]; original_pos=next(iter(mapped))
            for f in freq[(chrom,original_pos)]:
                if str(f['variant_y']) not in variants: continue
                if (ea,oa)==(f['ea_y'],f['oa_y']): matches.append((f['eaf_y'],f['source_y'],'exact_original_alleles'))
                elif (ea,oa)==(f['oa_y'],f['ea_y']): matches.append((1-f['eaf_y'],f['source_y'],'reversed_original_alleles'))
            if len(matches)==1 and pd.notna(matches[0][0]):
                row.update(eaf_y=matches[0][0],eaf_source=matches[0][1],eaf_attachment_status=matches[0][2],frequency_pos37=original_pos)
            else: row['eaf_attachment_status']='missing_or_ambiguous_exact_frequency'
        records.append(row)
    return pd.DataFrame(records)

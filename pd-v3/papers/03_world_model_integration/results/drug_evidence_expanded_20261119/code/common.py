import re, json, numpy as np, pandas as pd
O='/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results/drug_evidence_expanded_20261119'
A='/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/'; S6=A+'six_criteria_20261003/'
UP={'AGONIST','PARTIAL AGONIST','ACTIVATOR','POSITIVE ALLOSTERIC MODULATOR','POSITIVE MODULATOR','OPENER','EXOGENOUS PROTEIN','EXOGENOUS GENE','FULL AGONIST','BIASED AGONIST'}
DOWN={'INHIBITOR','ANTAGONIST','ALLOSTERIC ANTAGONIST','BLOCKER','NEGATIVE MODULATOR','NEGATIVE ALLOSTERIC MODULATOR','INVERSE AGONIST','GATING INHIBITOR','COMPLEMENT INHIBITOR','DISRUPTING AGENT'}
SILENCE={'ANTISENSE INHIBITOR','RNAI INHIBITOR','DEGRADER','GENE EDITING NEGATIVE MODULATOR'}
G_DOWN={'INHIBITION','ANTAGONIST','NEGATIVE','INVERSE AGONIST','PORE BLOCKER','COMPETITIVE','IRREVERSIBLE INHIBITION','FEEDBACK INHIBITION','NON-COMPETITIVE','VOLTAGE-DEPENDENT INHIBITION'}
G_UP={'FULL AGONIST','AGONIST','PARTIAL AGONIST','POSITIVE','ACTIVATION','POTENTIATION','BIASED AGONIST','IRREVERSIBLE AGONIST'}
def action_direction(a):
    """Returns (direction, is_silencing). Only the mechanism vocabulary is used, never the drug name."""
    a=str(a or '').upper().strip()
    if a.startswith('GTOP|'):
        _,t,x=(a.split('|')+['',''])[:3]
        if x in G_UP and t!='ANTAGONIST': return 'up',False
        if x in G_DOWN: return 'down',False
        if x=='PARTIAL AGONIST': return 'unknown',False
        if t in('INHIBITOR','ANTAGONIST','CHANNEL BLOCKER','GATING INHIBITOR') and x in('','BINDING','SLOWS INACTIVATION'): return ('down',False) if x!='SLOWS INACTIVATION' else ('unknown',False)
        if t in('AGONIST','ACTIVATOR') and x in('','BINDING'): return 'up',False
        return 'unknown',False
    if a.startswith('DGI|'):
        _,t,d=(a.split('|')+['',''])[:3]; ts=set(x for x in t.split(';') if x)
        real=ts-{'ANTIBODY','IMMUNOTHERAPY','BINDER','OTHER/UNKNOWN','VACCINE','MODULATOR'}
        if not real: return 'unknown',False
        if real<= {'INHIBITOR','ANTAGONIST','BLOCKER','NEGATIVE MODULATOR'}: return 'down',False
        if real<= {'AGONIST','ACTIVATOR','POSITIVE MODULATOR'}: return 'up',False
        return 'unknown',False
    if a in SILENCE: return 'down',True
    if a in UP: return 'up',False
    if a in DOWN: return 'down',False
    return 'unknown',False
DEPL=re.compile(r'antibody[- ]drug conjugate|immunoconjugate|radio(?:immuno|conjugate|label|therapeutic)|immunotoxin|chimeric antigen receptor|car[- ]?t\b|t[- ]cell engag|bispecific t|cytolytic|cytotoxicity|\badcc\b|\bcdc\b|deplet',re.I)
NEG=re.compile(r"non[- ]?(?:depleting|cytotoxic|cytolytic|lytic)|(?:without|lacks?|lacking|no|not|devoid of|reduced|minimal)\s[^.;]{0,40}?(?:deplet\w*|cytotox\w*|cytolytic|adcc|cdc)",re.I)
def depl_hit(t):
    return DEPL.search(NEG.sub(' ',str(t)))
def antigens():
    r=[x for x in json.load(open(A+'candidate_registry.json')) if not x['is_control']]
    return r
def norm_list(s): return [x for x in str(s).split(';') if x and x!='nan']

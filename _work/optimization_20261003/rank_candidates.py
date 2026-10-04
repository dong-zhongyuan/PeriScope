"""Join regenerated evidence; full table and explicitly gated priority lists."""

# Current authoritative pipeline dispatch (imports retain helper compatibility).
if __name__ == '__main__':
    from pathlib import Path as _RegistryPath
    import subprocess as _RegistrySubprocess, sys as _RegistrySys
    if _RegistryPath('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003/model_registry.json').exists():
        _RegistrySubprocess.run([_RegistrySys.executable, str(_RegistryPath(__file__).with_name('refresh_model_outputs.py'))], check=True)
        _RegistrySubprocess.run([_RegistrySys.executable, str(_RegistryPath(__file__).with_name('refresh_current_pipeline.py'))], check=True)
        raise SystemExit(0)

import json
from pathlib import Path
import numpy as np
import pandas as pd

O=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
RULES=dict(
 technical='unique gene mapping or defined protein complex, held-out donor-mean Spearman calculated within fine CITE-seq subtypes (equal subtype weights within each model blood state) >=0.20, at least two held-out donors with rho>0, nonzero measured dose IQR',
 reproducible='>=4/5 confirmation seeds in frozen program direction, >=0.75 adjacent dose steps monotone',
 specificity='gene-decoy q<=0.05 across all programs and protein-decoy percentile>=0.95; protein calibration is an effect-size criterion, not an adjusted significance claim',
 brain='Ma/Kamath competitive meta q<=0.05 across programs and disease directions; same relative enrichment sign and same actual donor-score disease direction in both cohorts',
 external='GSE157783 donor program contrast agrees with the actual Ma/Kamath donor-score disease direction; exact p and corrected q retained separately',
 tier_A='all above plus GWAS q<=0.05',tier_B='all above without requiring GWAS',
 tier_C='technically eligible, reproducible, gene-decoy q<=0.05 and brain criterion; protein-specificity or independent direction incomplete',
 ranking='tier, confirmation direction fraction, protein-decoy q, gene-decoy q, brain-meta q, held-out protein rho; no fitted composite weights',
 tier_C_status='follow-up candidates, not promoted first-choice targets')

def main():
    (O/'ranking_rules.json').write_text(json.dumps(RULES,indent=2))
    r=pd.read_csv(O/'randko_all_programs.csv');registry=json.loads((O/'candidate_registry.json').read_text())
    annotations={v['target']:v for v in registry}
    ac=pd.concat([pd.read_csv(p) for p in sorted((O/'evaluation').glob('protein_accuracy_seed*.csv'))],ignore_index=True)
    # Donors, not cells or seeds, are the units summarized for bridge support.
    ac=ac.groupby(['target','blood_state','donor'],as_index=False)[['rho','mse']].mean()
    bridge={}
    for (t,s),df in ac.groupby(['target','blood_state']):bridge[t,s]=dict(protein_rho=float(df.rho.mean()),protein_mse=float(df.mse.mean()),positive_protein_donors=int((df.rho>0).sum()),n_protein_donors=len(df))
    support=pd.read_csv(O/'all_candidate_input_support.csv').set_index(['target','blood_state'])
    brain=json.loads((O/'brain_spatial_competitive.json').read_text())
    ext=pd.read_csv(O/'independent_validation_stats.csv').set_index('combo')
    blood_assoc=pd.read_csv(O/'blood_disease_associations.csv')
    blood_assoc=blood_assoc[blood_assoc.cohort.astype(str)=='combined'].set_index(['target','blood_state','assay'])
    from program_effect_sizes import main as effect_sizes
    effect_sizes()
    effects=pd.read_csv(O/'program_effect_sizes.csv').set_index('combo')
    spatial=pd.read_csv(O/'spatial_localization.csv')
    spatial_by={c:d for c,d in spatial.groupby('combo')}
    assays=pd.read_csv(O/'all_target_assay_coverage.csv').set_index('target') if (O/'all_target_assay_coverage.csv').exists() else pd.DataFrame()
    drugs=pd.read_csv(O/'drug_gene_interactions.csv') if (O/'drug_gene_interactions.csv').exists() else pd.DataFrame()
    rows=[]
    for rec in r.to_dict('records'):
        combo=rec['combo'];p=rec['protein'];bs,br=rec['axis'].split('__x__');a=annotations[p]
        row=rec|dict(genes=';'.join(a['genes']),mapping=a['mapping'],blood_state=bs,brain_state=br)
        for assay in ['observed_RNA','inferred_ADT']:
            key=(p,bs,assay)
            if key in blood_assoc.index:
                ar=blood_assoc.loc[key];row['blood_'+assay+'_effect']=float(ar.effect_PD_minus_control_adjusted);row['blood_'+assay+'_q']=float(ar.q_BH)
        if p in assays.index:
            ar=assays.loc[p];row.update(Olink_gene_level_coverage=bool(ar.Olink_gene_level_coverage),Olink_covered_genes=ar.Olink_covered_genes)
        if len(drugs):
            dr=drugs[drugs.gene.isin(a['genes'])];row.update(n_drug_names=int(dr.drug.nunique()),n_approved_drug_names=int(dr.loc[dr.drug_approved==True,'drug'].nunique()))
        if combo in spatial_by:
            sr=spatial_by[combo];row.update(spatial_median_moran_I=float(sr.moran_I.median()),spatial_positive_FDR_sections=int(((sr.moran_I>0)&(sr.q_all_program_sections<=.05)).sum()),spatial_n_sections=len(sr),spatial_context_rho=float(sr.context_marker_rho.median()),spatial_context_rho_excluding_shared_genes=float(sr.context_rho_excluding_shared_genes.median()))
        row.update(bridge.get((p,bs),{}));row['dose_iqr']=float(support.loc[(p,bs),'dose_iqr'])
        sid='gsea::'+combo
        candidates=[(brain['meta_stouffer'].get(d,{}).get(sid,{}).get(f'q_{d}_meta',1),d) for d in ['p_up','p_down']]
        q,d=min(candidates);sign=1 if d=='p_up' else -1
        ma=brain['ma_sets'].get(sid,{});ka=brain['kamath_sets'].get(sid,{})
        ma_effect=ma.get('stat',np.nan)-ma.get('null_mean',np.nan);ka_effect=ka.get('stat',np.nan)-ka.get('null_mean',np.nan)
        ma_actual=ma.get('donor_score_effect_PD_minus_control',np.nan);ka_actual=ka.get('donor_score_effect_PD_minus_control',np.nan)
        actual_sign=float(np.sign(ka_actual)) if np.isfinite(ma_actual) and np.isfinite(ka_actual) and ma_actual*ka_actual>0 else 0
        row.update(brain_meta_q=q,competitive_direction='higher_than_matched_background' if sign==1 else 'lower_than_matched_background',disease_direction='up' if actual_sign>0 else 'down' if actual_sign<0 else 'unresolved',ma_donor_score_effect=ma_actual,kamath_donor_score_effect=ka_actual,ma_relative_effect=ma_effect,kamath_relative_effect=ka_effect,
                   ma_q=ma.get(f'q_{d}_ma',1),kamath_q=ka.get(f'q_{d}_kam',1))
        er=ext.loc[combo].to_dict() if combo in ext.index else {}
        row.update(independent_effect=er.get('effect_PD_minus_control',np.nan),independent_p=er.get('p_exact_two_sided',np.nan),independent_q=er.get('q_BH_all_programs',np.nan))
        mapping_pass=a['mapping'] in ['unique','defined_complex']
        protein_pass=row.get('protein_rho',-1)>=.2 and row.get('positive_protein_donors',0)>=2
        dose_pass=row['dose_iqr']>1e-6
        technical=mapping_pass and protein_pass and dose_pass
        row.update(mapping_pass=mapping_pass,protein_prediction_pass=protein_pass,measured_dose_pass=dose_pass)
        reproducible=rec['confirmation_positive_fraction']>=.8 and rec['confirmation_monotonic_fraction']>=.75
        specificity=rec['q_gene_decoy']<=.05 and rec['protein_percentile']>=.95
        brain_pass=q<=.05 and sign*ma_effect>0 and sign*ka_effect>0 and actual_sign!=0
        external=actual_sign*row['independent_effect']>0
        tier='D'
        if technical and reproducible and rec['q_gene_decoy']<=.05 and brain_pass:tier='C'
        if technical and reproducible and specificity and brain_pass and external:tier='B'
        row.update(tier=tier,technical_pass=technical,reproducible_pass=reproducible,specificity_pass=specificity,brain_pass=brain_pass,independent_direction_pass=external)
        row['unmet_nomination_criteria']=';'.join(name for name,passed in [('gene_or_complex_identity',mapping_pass),('held_out_protein_prediction',protein_pass),('measured_dose_range',dose_pass),('confirmation_seeds',reproducible),('specificity',specificity),('brain_disease_support',brain_pass),('independent_direction',external)] if not passed)
        if combo in effects.index:
            row.update(effects.loc[combo].to_dict())
        else:row['predicted_restorative_protein_change']='unresolved'
        rows.append(row)
    df=pd.DataFrame(rows) if rows else pd.DataFrame(columns=['combo','protein','genes','tier'])
    if len(df):df=df.sort_values(['tier','confirmation_positive_fraction','q_protein_decoy','q_gene_decoy','brain_meta_q','protein_rho'],ascending=[True,False,True,True,True,False])
    df.to_csv(O/'all_program_evidence.csv',index=False)
    first=df[df.tier.isin(['A','B'])].drop_duplicates('protein') if len(df) else df
    first.to_csv(O/'new_priority_targets.csv',index=False)
    (first.drop_duplicates('genes') if len(first) else first).to_csv(O/'new_priority_genes.csv',index=False)
    best=df.drop_duplicates('protein') if len(df) else df
    full=[]
    for a in registry:
        if a['is_control']:continue
        hit=best[best.protein==a['target']] if len(best) else best
        rec=hit.iloc[0].to_dict() if len(hit) else dict(protein=a['target'],tier='unannotated_response',reason='no discovery-defined enriched program',genes=';'.join(a['genes']),mapping=a['mapping'])
        full.append(rec)
    pd.DataFrame(full).to_csv(O/'all_targets.csv',index=False)
    summary=dict(n_biological_antigens=len(full),n_programs=len(df),n_priority_targets=len(first),tier_counts=df.tier.value_counts().to_dict() if len(df) else {},
                 priority_targets=first.protein.tolist() if len(first) else [],complete=True)
    (O/'rescreen_summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary),flush=True)

if __name__=='__main__':
    main()
    from validate_outputs import main as validate
    validate()

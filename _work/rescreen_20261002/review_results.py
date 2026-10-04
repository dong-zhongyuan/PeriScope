"""Summarize complete rerun evidence without changing the frozen ranking."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd

O=Path('/public/home/mengxl/dzy/pd_product_assets/results/rescreen_20261002')

def main():
    qa=json.loads((O/'final_quality_checks.json').read_text());assert qa['status']=='passed'
    r=pd.read_csv(O/'all_program_evidence.csv');alltargets=pd.read_csv(O/'all_targets.csv');priority=pd.read_csv(O/'new_priority_targets.csv')
    programs=json.loads((O/'programs.json').read_text());ab=json.loads((O/'ablation_comparison_definition.json').read_text())
    assert ab['programs_sha256']==hashlib.sha256((O/'programs.json').read_bytes()).hexdigest()
    condition=np.ones(len(r),bool);flow=[]
    for label,key in [('全部发现程序',None),('抗原身份', 'mapping_pass'),('留出供体蛋白预测','protein_prediction_pass'),('测量剂量范围','measured_dose_pass'),('确认种子稳定性','reproducible_pass'),('randKO选择性','specificity_pass'),('脑疾病支持','brain_pass'),('独立队列方向','independent_direction_pass')]:
        if key is not None:condition &= r[key].fillna(False).to_numpy(dtype=bool)
        flow.append(dict(step=label,n_programs=int(condition.sum()),n_antigens=int(r.loc[condition,'protein'].nunique())))
    pd.DataFrame(flow).to_csv(O/'nomination_flow.csv',index=False)
    target_flow=[]
    for target,d in r.groupby('protein'):
        keep=np.ones(len(d),bool)
        row=dict(protein=target,n_programs=len(d))
        for key in ['technical_pass','reproducible_pass','specificity_pass','brain_pass','independent_direction_pass']:
            keep &= d[key].fillna(False).to_numpy(dtype=bool)
            row['remaining_after_'+key]=int(keep.sum())
        target_flow.append(row)
    pd.DataFrame(target_flow).to_csv(O/'target_nomination_flow.csv',index=False)
    folds=json.loads((O/'program_folds.json').read_text());fold_rows=[]
    for name,parts in folds.items():
        fa=set(parts.get('A',[]));fb=set(parts.get('B',[]))
        fold_rows.append(dict(combo=name,n_fold_A=len(fa),n_fold_B=len(fb),n_shared=len(fa&fb),n_union=len(fa|fb),jaccard=len(fa&fb)/max(len(fa|fb),1),present_both_folds=bool(fa and fb)))
    fold_table=pd.DataFrame(fold_rows);fold_table.to_csv(O/'discovery_program_reproducibility.csv',index=False)
    fold_lookup=fold_table.set_index('combo')
    gene=pd.read_csv(O/'target_gene_expression_by_donor.csv')
    # Retain raw biological context alongside the statistical nomination, not as a tuned new gate.
    context=[]
    selected=pd.concat([priority,alltargets[alltargets.protein.isin(['CD22','CD35'])]],ignore_index=True).drop_duplicates('protein')
    for v in selected.to_dict('records'):
        bs=v.get('blood_state');genes=str(v.get('genes','')).split(';');d=gene[(gene.gene.isin(genes))&(gene.blood_state==bs)]
        row=dict(protein=v['protein'],genes=v.get('genes'),blood_state=bs,tier=v.get('tier'),combo=v.get('combo'),unmet_nomination_criteria=v.get('unmet_nomination_criteria',''))
        row.update(n_mapped_genes=len(genes),n_mapped_genes_in_model_RNA=int(d.gene.nunique()),mapped_genes_in_model_RNA=';'.join(sorted(set(d.gene))))
        if v.get('combo') in fold_lookup.index:row.update(fold_lookup.loc[v['combo']].to_dict())
        if len(d):
            row.update(observed_blood_gene_detection_median=float(d.fraction_nonzero.median()),observed_blood_gene_detection_max=float(d.fraction_nonzero.max()),observed_blood_mean_log1p_CP10k=float(d.mean_log1p_CP10k.mean()),n_observed_blood_donors=d.donor.nunique())
        context.append(row)
    pd.DataFrame(context).to_csv(O/'priority_and_previous_target_context.csv',index=False)
    ext=pd.read_csv(O/'independent_validation_stats.csv');spatial=pd.read_csv(O/'spatial_localization.csv');tissue=pd.read_csv(O/'tissue_axis_associations.csv')
    cov=pd.read_csv(O/'blood_disease_covariate_sensitivity.csv');a=pd.read_csv(O/'ablation_program_response_comparison.csv')
    summary=dict(status='complete',n_antigens=qa['n_biological_antigens'],n_programs=len(programs),
        n_priority_antigens=len(priority),n_priority_gene_groups=len(priority.drop_duplicates('genes')),priority_targets=priority.protein.tolist(),
        nominated_program_tiers=r.tier.value_counts().to_dict(),nomination_flow=flow,
        independent_n_tested=len(ext),independent_FDR05=int((ext.q_BH_all_programs<=.05).sum()),
        discovery_both_fold_programs=int(fold_table.present_both_folds.sum()),discovery_median_gene_jaccard=float(fold_table.jaccard.median()),
        spatial_n_sections=int(spatial['sample'].nunique()),spatial_n_program_sections=len(spatial),
        tissue_associations=len(tissue),blood_covariate_sensitivity_rows=len(cov),
        ablations={k:dict(n_programs=len(d),median_response_norm_retention=float(d.response_norm_retention.median()),median_response_cosine=float(d.response_cosine.median())) for k,d in a.groupby('control')},
        auxiliary_checks=dict(ablation_program_hash=True,drug_registry_hash=json.loads((O/'drug_annotation_source.json').read_text()).get('registry_sha256')==hashlib.sha256((O/'candidate_registry.json').read_bytes()).hexdigest()))
    assert all(summary['auxiliary_checks'].values())
    (O/'review_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)

if __name__=='__main__':main()

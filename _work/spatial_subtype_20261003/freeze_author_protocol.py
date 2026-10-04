from pathlib import Path
import json
P=Path('/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003')
O=P/'spatial_subtype_author_20261003';O.mkdir(exist_ok=True)
assert not (O/'samples').exists()
p=json.loads((P/'spatial_subtype_20261003/protocol.json').read_text())
p.update(version='20261003_spatial_subtype_author_v6',reference='Ma et al author annotated substantia nigra RNA counts,196693 nuclei; GEO identity crosswalk; exclude all overlapping spatial donors and unresolved OJ5 identity; donor-balanced reference',
 states=['astro_fibrous','astro_protoplasmic','microglia_homeostatic','microglia_mhc2'],
 recipient_mapping={'astro':['astro_fibrous','astro_protoplasmic'],'microglia_homeostatic':['microglia_homeostatic'],'microglia_mhc2':['microglia_mhc2']},
 resolution='Author annotated Fibrous/Protoplasmic astrocytes and Baseline/Activated microglia; Monocyte_Derived BAM separate; model astro parent programs evaluated in both spatial leaves, without claiming the model predicts leaf identity',
 reference_QC_holdout=['OJ11','OJ13','OJ9','OJ22'],
 deconvolution_genes='40 reference enriched markers per author subtype plus pre-existing biological markers; exclude sex/MT/ribosomal genes; evaluation genes disjoint',
 amendment='Kamath-derived microglia references were not identifiable in heldout mixtures; replace with original spatial-study author subtype annotations and donor-disjoint reference before any candidate outcome inspection',
 leaf_aggregation='BH over all program x applicable spatial leaf x both directions; use q in pre-existing matched-brain direction AND same actual donor effect sign; any supported leaf passes, no support plus unresolved leaf remains unresolved; all tested nonsupport means not_supported',
 candidate_outcomes_used_in_reference_selection=False)
(O/'protocol.json').write_text(json.dumps(p,indent=2))

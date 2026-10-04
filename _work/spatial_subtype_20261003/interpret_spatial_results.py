from pathlib import Path
import json,os
import pandas as pd
P=Path(os.environ.get('PERISCOPE_RESULT_ROOT','/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003'))
S=Path(os.environ.get('SPATIAL_OUTPUT_DIR',str(P/'spatial_subtype_author_20261003')));O=P/'seven_criteria_20261003'
e=pd.read_csv(O/'all_6849_programs_seven_criteria.csv');summary=json.loads((O/'completed.json').read_text())
pos=e[e.C7.eq('pass')];bio=e[e.five_criteria_pass&e.C7.eq('pass')];a=e[e.seven_criteria_pass]
lines=['空间亚型验证：结果解读与采用建议','',
'已完成原始空间计数的亚型分解、逐供者亚型表达估计和213个抗原/6849条程序的重筛。星形胶质拆为Fibrous与Protoplasmic，小胶质拆为作者Baseline与Activated（MHC-II富集）状态；BAM独立建模。',
'现有模型的astro仍是父级类型，本次空间分析将其程序定位到实际作者亚型，未把模型输出虚称为更细亚型预测。','',
'本次得到的积极结果：',
f'1. 空间本身支持{pos.protein.nunique()}个抗原、{len(pos)}条程序。',
f'2. 原前五项生物学标准+空间支持：{bio.protein.nunique()}个抗原、{len(bio)}条程序，包括'+ '、'.join(sorted(bio.protein.unique()))+'。',
f'3. 再满足原第六项直接临床药物条件，即七项同时通过：{a.protein.nunique()}个抗原、{len(a)}条程序。',
'4. 七项通过者CD80的两条程序，在旧混合组织结果中与脑细胞方向相反；在Fibrous星形胶质内与原脑细胞方向一致。说明用混合组织方向强制否决细胞状态程序确实可能错失候选。','',
'CD80逐程序证据：']
for r in a.itertuples():
 lines.append(f'{r.blood_state} → {r.spatial_subtype}: {r.n_PD} PD/{r.n_control}对照；覆盖{r.n_used}/{r.n_program_genes}个程序基因；空间程序得分PD减对照={r.effect:.4f}；表达匹配竞争检验q={r.spatial_q:.4f}；精确供者标签置换p={r.spatial_donor_permutation_p:.4f}。')
lines += ['',
'该信号是相对表达匹配背景的程序富集，并伴有与脑细胞一致的疾病方向；整体供者得分差异目前未达显著。药物记录表示靶点已有直接临床药物关联，不表示这些药物在PD中的效果已经得到本项目验证。','',
'采用方式：',
'保留原六项作为全局候选入口，把七项同时通过者作为新增的空间支持层。尚未可靠区分的亚型和基因覆盖不足的程序保留待判定，不能等同阴性。',
'技术覆盖：Fibrous星形胶质5 PD/5对照、MHC-II小胶质3 PD/4对照可用。Protoplasmic星形胶质的独立恢复较弱；稳态小胶质还缺足够可估供者。',
'原六项356条程序的去向：2条支持、153条已检验未支持、201条暂不可判定。同一抗原可有多条程序，因此不同状态下的抗原数不能简单相加。',
'部署检查中修正了父级/子型的覆盖判断：父级下任一可靠亚型的阳性证据可以支持一条程序；若另一个相关亚型仍不可判定，不能据此对整个父级作阴性淘汰。修正前后全部程序统计量、C7状态与七项名单的文件哈希完全一致，记录见adoption_scope_clarification.json。','',
'CD22：原六项通过的cDC→MHC-II小胶质程序本轮竞争q约0.723；cDC→astro程序基因覆盖不足。',
'CD35：本轮未进入七项支持层，既有直接药物条件仍未建立；所有程序结果已保留。','',
'文件：',
'七项通过候选汇总.csv / A7_all_programs.csv：严格七项支持层。',
'five_biological_plus_spatial_without_drug_requirement.csv：前五项+空间的10个抗原。',
'six_to_seven_changes.csv：原六项每条程序的去向。',
'A6_pending_spatial_resolution.csv：原六项中201条暂不可判定程序。',
'all_leaf_program_spatial_tests.csv（空间目录内）：全部9308个亚型检验，无选择性删除。',
'本轮不修改论文、主图或原六项门槛。']
(O/'结果解读与采用建议.txt').write_text('\n'.join(lines)+'\n')
print('INTERPRETATION_SAVED')

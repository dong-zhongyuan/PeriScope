"""gen_layouts.py - single source of the five page layouts. Writes layout_fig1.json ... layout_fig5.json and LAYOUT.md.
Block = rows x cols of the 4 x 5 unit grid (u = 41.25 mm, g = 3 mm); row, col = 1-based top-left cell."""
import json
import os
import sys

OUT = sys.argv[1] if len(sys.argv) > 1 else '.'
U, G = 41.25, 3.0
ABBR = {
    'R03': '/public/home/mengxl/dzy/pd-v3/papers/03_world_model_integration/results',
    'R01': '/public/home/mengxl/dzy/pd-v3/papers/01_mr_imaging_scan/results',
    'OPT': '/public/home/mengxl/dzy/pd_product_assets/results/optimization_20261003',
    'FRD': '/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design/manuscript_review/figure_redesign_20261004/data',
    'DES': '/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/project/pd_product_figure_design',
    'NMI': '/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/delivery/journal_versions/nature_machine_intelligence',
    'PIG': '/public/home/mengxl/dzy/pd_product_assets/results/pig_external_validation_20261003',
    'MOUSE': '/public/home/mengxl/dzy/pd_product/manuscript_revision_20261017/diag/periscope_panel_correlation_20261110/payload/pd_product_figure_design/manuscript_review/periscope_panel_correlation_20261110/mouse_followup',
    'HH': '/public/home/mengxl/dzy/pd-v3/papers/figures_v2_20261119/_system/reuse_hh (code copied from the Mac hh/papers/01b_breast_us_algo/figures_r)',
}
READY, BLANK = 'ready', 'blank-pending-ablation'

FIGS = [
    (1, '模型构建', [
        ('a', '1a_architecture', '模型架构：血液 RNA 编码器（锁定 64 维血液状态）→ 228 通道蛋白头 → 耦合势能 → 脑 RNA 解码器，Sinkhorn 监督；框内印 813,363 个观测、228 个 ADT、4,741 个基因；右端并入留出 CITE-seq 评估支路（P6–P8 只输入 RNA、模型冻结、同细胞 ADT 作参照）', '1x3', 1, 1,
         'NMI figures/ 与 pptx/ 的原 1A、2A 图标；DES/workflow_assets/blender/Figure1A_inventory.csv', '原 1A 图标；hh 图 1a 数量入框（HH/Fig1/A.R）；hh 图 5a 编号步骤', READY, '1a + 1g'),
        ('b', '1b_resources', '资源构成嵌套三环：内环观测类型、中环队列、外环队列内诊断；洞内印总数 813,363；外圈 10% 刻度；环形图标作图例，各类印占比', '2x1', 1, 4,
         'DES/full_template_reuse/inputs/all_cohort_composition.csv；DES/template_semantic_review/data/input_1A.csv', '#01 aivc NAR 1B；dring()、ring_ticks()、ring_icon()', READY, '1b'),
        ('c', '1c_passports', '六个转录组资源加 CITE-seq 桥的"护照"：每行一个资源，对齐轨道依次为组织图标、对照与 PD 供体单位方块、对数细胞数、训练或留出用途、每供体细胞数蜂群（并入原 2C 留出细胞数）', '1x2', 2, 1,
         'DES/full_template_reuse/inputs/all_cohort_composition.csv；FRD/protein_coverage.csv、protein_bridge.csv', '#02 hh 图 1c（HH/Fig1/C.R）；waffle()', READY, '1f'),
        ('d', '1d_blood_mosaic', '6 个血液供体（Control、Early PD、Late PD）的细胞类型构成马赛克：条宽 = 供体细胞数，段高 = 细胞类型比例', '1x1', 2, 3,
         'DES/full_template_reuse/inputs/all_cohort_composition.csv（GSE223138）', 'hh S2d 马赛克（HH/S2/D.R）', READY, '1d'),
        ('e', '1e_brain_umap', '脑 RNA UMAP 两张：灰度 hexbin = 全部训练细胞密度，上叠彩色留出细胞；一张按诊断着色，一张按细胞类型着色（类型名直接标在簇中心）', '1x2', 3, 1,
         'DES/panel_grid/objects/1C.rds 及 umap_grid.R 的输入（原 1C）', 'hh 图 1e 密度底加叠点（HH/Fig1/E.R，hex_bins()）', READY, '1c'),
        ('f', '1f_annotation_alluvial', '已发表血细胞注释到建模细胞状态（cDC、cMono、nMono、pDC 等）的冲积图，两端印细胞数；右侧并排标记基因点热图（色 = 均值，大小 = 检出率；原 1E 并入）', '1x2', 3, 3,
         'DES/panel_grid/objects/1E.rds、1F.rds 的源数据（原 1E、1F）', '#03 CREED NC 1f、1e（creed/figures_nc/r_Figure1/F.R、E.R）；dribbon()', READY, '1e'),
        ('g', '1g_training_loss', '蛋白头训练损失：10 个种子（42–51）细线 + 均值粗线，停止点画圆点，线端印末端损失', '1x1', 4, 1,
         'FRD/training.csv', '#04 HepaWorld S1a', READY, '1h'),
        ('h', '1h_protein_mse', '留出供体 × 亚型的蛋白 MSE 逐单位斜率线：PeriScope 对 Ridge、State mean、Subtype mean；"总体"与"亚型平衡"两个分面，均值粗线', '1x2', 4, 2,
         'FRD/protein_accuracy.csv、protein_accuracy_model_state_seed42.csv（原 2D 源表）', '#10 hh 图 6e（HH/Fig6/E.R）', READY, '1i'),
        ('i', '1i_subtype_spearman', '非对照蛋白的亚型内 Spearman：5 个亚型 × 3 个供体共 15 行，每行一个蛋白一个点（10 种子平均）的蜂群 + 中位短线 + 右侧中位数列', '2x1', 4, 4,
         'FRD/protein_accuracy.csv', '#14 CREED paper2 图 2c 行设计', READY, '1k'),
        ('j', '1j_final_target_adt', '最终靶点 CD71/TFRC、CD119/IFNGR1 的预测 ADT 对实测 ADT：2 个靶点 × P6、P7、P8 共 6 格 hexbin，每格印 n 与 rho，画拟合线与恒等线', '1x2', 5, 1,
         'FRD/protein_cell_pairs.csv（靶点换为 CD71、CD119）', '#05 hh S1b（HH/S1/B.R + hex_bins()）', READY, '1j'),
        ('k', '1k_antigen_consistency', '逐抗原一致性：228 个抗原按 Spearman 排序，带 10 种子区间；中位 0.63 虚线；最终靶点标名；右侧边缘密度', '1x1', 5, 3,
         'FRD/protein_accuracy.csv 派生的逐抗原 Spearman（源表请主会话确认）', '#20 hh 图 4g 排序分数条（HH/Fig4/G.R）', READY, '1l'),
    ]),
    (2, '模型能力', [
        ('a', '2a_design', '评估设计：4 个血液状态 × 2 个脑状态 = 8 条路线的二部图；测试层级图标（人类疾病供体、外部队列、猪 LPS、小鼠 LPS）；10 个种子', '1x1', 1, 1,
         'NMI panels/3a.svg（原 3A）', '原 3A 图标', READY, '2a'),
        ('b', '2b_leaderboard', '主视觉。人类疾病方向 Spearman 排行榜瓦片：方法 × 8 条路线 + 均值列；填色 = 种子均值（PeriScope 10 个种子，比较方法按各自种子数）；瓦片内逐种子小点（不用小柱，比较方法只有 3 个种子）；气泡 = 胜过 CellOT 的种子数；低于 CellOT 的瓦片加斜线；右列相对提升（比 CellOT 高 61%，比 OT-Ridge 高 176%）', '2x2', 1, 2,
         'R03/corrected_benchmark_20261119/TABLE_benchmark.csv', '#06 aivc NAR 3B；备选 #07 CREED NC 2c', READY, '2b'),
        ('c', '2c_pearson_direction', 'Pearson r（填色）与方向一致率（点大小）点热图：方法 × 路线 + 均值列，方法行与 b 逐行对齐', '2x1', 1, 4,
         'R03/corrected_benchmark_20261119/TABLE_benchmark.csv', '#03 右半点热图；#07', READY, '2c'),
        ('d', '2d_seed_stability', '种子稳定性：8 条路线各一行，10 个种子（42–51）蜂群点，最佳比较方法均值画刻度，右列"胜过最佳比较方法的种子数 k/10"；底行为逐个去掉 locked 测试供体的范围', '1x1', 2, 1,
         'R03/corrected_benchmark_20261119/TABLE_ablation_corrected_by_seed.csv（只取完整模型行）、TABLE_locked_leave_one_test_donor.csv', '#14 行设计；YeastBridge 图 6 记分卡 Met 列', READY, '2h'),
        ('e', '2e_external_cohorts', '外部人类队列三张卡（GSE157783、GSE184950、Ma 2025）：卡头色块印相对 CellOT 的倍数（2.8×、2.1×，Ma 如实画出略低）；卡内每个方法一行，donor bootstrap 密度带 + 点估计；灰带为置换零分布', '1x3', 3, 1,
         'R03/corrected_benchmark_20261119/TABLE_benchmark_key_summary.csv 及同目录逐队列 bootstrap、置换表', '#08 hh 图 3b（HH/Fig3/B.R + draw_strip()、ci_rail()）；备选 #09', READY, '2d'),
        ('f', '2f_readout_optimisation', '读出优化：原始读出 → 10 种子集成 → 方差匹配读出三步；纵轴为相对原始读出的距离（%），MMD 与 energy × Astro、MG-AP 四条线，原始步画 10 个单种子点；终点印降幅（MMD 53–54%，energy 76%）。只用模型结果', '1x1', 3, 4,
         'R03/metric_optimization_20261119/distribution_adt（README 主表，仅模型行）', '#13 hh 图 6b 的分解思路，压缩为相对值阶梯', READY, '2i'),
        ('g', '2g_structural_ablation', '结构消融（完整模型、随机配对、去耦合）逐路线斜率线', '1x2', 4, 1,
         'R03/manuscript_readout_ablation_20261119（新消融评估运行中）', '#10 hh 图 6e', BLANK, '2e'),
        ('h', '2h_direction_ablation', '程序响应方向消融：消融臂 × 路线影响矩阵', '1x1', 4, 3,
         '新消融评估（运行中）', '#11 aivc NAR 6B', BLANK, '2f'),
        ('i', '2i_distribution_ablation', '分布指标消融（MMD、energy、均值 MSE）', '1x1', 4, 4,
         'R03/corrected_benchmark_20261119/TABLE_ablation_by_readout.csv 等（待新评估）', '#12 hh 图 2b', BLANK, '2g'),
        ('j', '2j_pig_lps', '猪 LPS 前额叶响应：左为 3,730 个直系同源基因预测对实测 hexbin（PeriScope）；右为方法森林，每行 1,000 次分组配对动物 bootstrap 密度带；最佳比较方法的 2.83×', '1x2', 5, 1,
         'PIG/（pig_summary、基因与 bootstrap 表）', '#15 hh 图 3c（HH/Fig3/C.R）；#05', READY, '2j'),
        ('k', '2k_mouse_lps', '小鼠小胶质 LPS 两个留出研究（Boles、Shemer）：左为 3,381 个基因预测对实测 hexbin 按研究分面；右为每只小鼠零响应参照到 PeriScope 的 MSE 配对线 + 研究均值（降 56.77%、51.55%）', '1x2', 5, 3,
         'MOUSE/', '#05；hh S5d 分位点图（HH/S5/D.R）', READY, '2k'),
    ]),
    (3, '模型读出', [
        ('a', '3a_readout_flow', '读出流程：五个剂量查询 → 发现种子 42–46 与确认种子 47–51 → 冻结程序（6,849 个）→ 匹配脑状态与空间检验', '1x1', 1, 1,
         'NMI 原 4A、5A 图标', '原 4A + 5A 图标', READY, '3a'),
        ('b', '3b_example_programs', '五个示例程序卡：上部剂量—响应曲线（种子细线，发现均值与确认均值粗线，卡头印斜率与方向）；下部该程序主要基因响应点（色 = 确认斜率 × 10³，大小 = 正斜率种子比例）', '1x3', 1, 2,
         'FRD/five_examples.csv、program_dose_curves.csv、gene_responses.csv', '#16 hh 图 4d（HH/Fig4/D.R）；#03 点热图', READY, '3b + 3c'),
        ('c', '3c_discovery_confirmation', '发现斜率对确认斜率：4,741 个基因 hexbin，极端基因标名', '1x1', 2, 1,
         'FRD/gene_responses.csv', '#05', READY, '3d'),
        ('d', '3d_gene_multiplicity', '每个基因出现在多少个程序（2,503 个基因）：对数横轴点直方图（茎加点）；307 个只在 1 个程序、64 个在 1,000 个以上', '1x1', 2, 2,
         'FRD/gene_reuse.csv', 'hh S2a 点直方图（HH/S2/A.R）', READY, '3j'),
        ('e', '3e_enrichment', 'CD22 负向查询成员径向基因星图 + CD40 正向查询 GO 点图（位置 = 重叠基因数，色 = −log10 q）', '1x2', 2, 3,
         'FRD/discovery_go.csv、cd40_functional_go.csv', '原 4E、4F', READY, '3e'),
        ('f', '3f_priority_programs', '最终靶点程序（CD71 4 个、CD119 6 个，共 10 行）× 10 个种子的带符号斜率发散热条（竖线分开发现与确认种子，右列均值）+ 表达匹配零分布小山脊与实测竖线 + 两种匹配方式的超额响应哑铃 + BH q 列；行序与 g 共用', '1x2', 3, 1,
         'OPT/program_seed_scores.csv、all_program_evidence.csv；R03/target_funnel_final_20261119/final_targets.csv（programs_gate1to5）', '#17 hh 图 6f（HH/Fig6/F.R + div_col()）；#19 对齐列', READY, '3f + 3g'),
        ('g', '3g_program_brain_support', '与 f 同行序：留一种子后程序是否保留（每行 5 格，q ≤ 0.05 实心）+ 匹配脑状态下 PD 减对照的供体效应（Kamath、GSE157783、Ma 2025 三点灰线相连，画零线）', '1x2', 3, 3,
         'FRD/matched_donor_support.csv、priority_seeds.csv；OPT/cell_state_sensitivity/', '#26 waffle()（HH/Fig5/C.R）；#19 对齐列', READY, '3h + 3i'),
        ('h', '3h_spatial_maps', '纤维型星形胶质与 MG-AP 的 RNA 分数空间图：2 行 × 10 张切片，切片按条件排列并加条件条，共用 0–1 平方根色标，缺失为灰；去卷积 rho 0.60 作数据标签', '1x4', 4, 1,
         'FRD/spatial_fractions.csv、spatial_quality.csv', '原 5G、5H；pd-v3 02 空间图', READY, '3m'),
        ('i', '3i_spatial_program_scores', '空间支持程序在纤维型星形胶质中的供体分数：程序 × 切片热图 + 条件条 + 证据状态列', '1x1', 5, 1,
         'FRD/spatial_heat_exact_matrix.csv、spatial_heat_annotation.csv、spatial_candidates.csv', '原 5C；hh S9a 热表（HH/S9/A.R）', READY, '3k'),
        ('j', '3j_donor_scores', '代表程序的供体分数（GSE157783，6 个对照、5 个 PD）：单位点 + 均值 + 95% 区间，三个分面，不用小提琴', '1x1', 5, 2,
         'FRD/matched_independent_donor_scores.csv', 'AgingAxis 图 6 蜂群', READY, '3l'),
        ('k', '3k_model_nomination', '模型侧初步提名单位账本：每行一个门（C1 身份与来源细胞、C2 剂量、C3 特异模型响应、C4 匹配细胞状态 PD 脑支持、C5 表面可及性），每格一个抗原（213 个），颜色表示通过或在此淘汰，最终靶点着色；右列印剩余数与占比（213 → 170 → 170 → 163 → 148 → 148）', '1x2', 5, 3,
         'R03/target_funnel_final_20261119/antigen_gate_matrix.csv、funnel_counts.csv', '#26 hh 图 5c（waffle()）；#20', READY, '3n'),
    ]),
    (4, 'MR 和共定位', [
        ('a', '4a_design', '设计：UKB-PPP 与 AGES 血浆 cis-pQTL（2,941 个 assay）→ 3,963 项影像表型 → 单工具 MR → 单信号 coloc 与 SuSiE-coloc → 反向 MR → BIG40 复现、cis-pQTL 复现，每步印数量', '1x1', 1, 1,
         'R01/full_atlas_20261116、R01/all_target_delivery_20261119', 'hh 图 5a 编号步骤；原 6A 图标', READY, '4a'),
        ('b', '4b_scale', '规模对比哑铃梯：assay、影像表型、MR 检验、区域、coloc 区域、反向 MR 六行，对数轴由旧值连到新值，右列印倍数（影像表型 248×，MR 检验约 20,000×）', '1x1', 1, 2,
         'R01/legacy_comparison.csv、R01/full_atlas_20261116', 'hh S2b 哑铃加数字列（HH/S2/B.R）', READY, '4b'),
        ('c', '4c_atlas', '全图谱：每个 assay 的最小 P 按 cis 工具所在染色体位置排列，UKB-PPP 向上、AGES 向下镜像', '1x2', 1, 3,
         'R01/full_atlas_20261116', '#21 figures_atlas_r（render.R + A.R）', READY, '4c'),
        ('d', '4d_category_heatmap', '类别效应谱：排名靠前的 assay × 影像类别最佳 z 发散热图，格内画 z 范围刻度，缺失画叉；CD71/TFRC、CD119/IFNGR1 所在行加框', '2x2', 2, 1,
         'R01/full_atlas_20261116', '#21（B.R）；R01/../figures_r/B.R', READY, '4d'),
        ('e', '4e_big40', 'BIG40 复现：左为发现 z 对复现 z hexbin（区域工具、MR 命中、匹配零分布三类）；右为同向率阶梯（最严一组 99.8%、全部命中 94.1%、同 z 带零分布，Wilson 区间）', '1x2', 2, 3,
         'R01/big40_repro_20261119/repro_summary.csv 及配对表', '#23 全靶点总览 g；#05', READY, '4e'),
        ('f', '4f_coloc_landscape', '单信号 coloc 全景（1,375 个区域）：上为 H4 沿基因组分布，下为两个窗口的后验条码热图', '1x2', 3, 3,
         'R01/all_target_coloc_20261118', '#22 全靶点总览 a、b（R01/all_target_overview_figure_20261119/scripts/render_overview.R）', READY, '4f'),
        ('g', '4g_susie', 'SuSiE-coloc：单信号 H4 对 SuSiE 最大 H4 的 hexbin（AGES、UKB-PPP 两个平台分面）+ 可信集是否包含 MR 工具的条码列', '1x2', 4, 1,
         'R01/all_target_susie_coloc_20261119/susie_coloc_pairs.csv.gz', '#05 hh S1b 加阈值格', READY, '4g'),
        ('h', '4h_reverse_mr', '反向 MR：左为四种策略下的配对估计散点；右为状态 × 染色体覆盖表，格内印数', '1x2', 4, 3,
         'R01/all_target_reverse_20261118；R01/all_target_delivery_20261119/reverse_status_table.csv', '#23 全靶点总览 f、h', READY, '4h'),
        ('i', '4i_final_target_forest', '最终靶点 CD71/TFRC、CD119/IFNGR1 × 关键影像表型的 MR 森林：左右侧实心与空心，模态着色；右侧 q、F 统计量、BIG40 同向三列', '1x2', 5, 1,
         'R03/target_funnel_final_20261119/final_target_mr_coloc_pairs.csv；R01/full_atlas_20261116', '#25 HepaWorld 4e–h 带统计列森林；#15', READY, '4i'),
        ('j', '4j_pqtl_replication', 'cis-pQTL 复现：全部 MR 工具在 UKB-PPP 与 AGES 中的效应散点，最终靶点工具加哑铃区间', '1x1', 5, 3,
         'R01/exposure_replication_20261116', '#25；#05', READY, '4j'),
        ('k', '4k_regional', '最终靶点区域关联：蛋白与影像两条 LocusZoom 式轨道，标出 SuSiE 可信集', '1x1', 5, 4,
         'R01/all_target_coloc_20261118、R01/all_target_susie_coloc_20261119 的区域表', '#24 pd-v3 plot_public_cd40.R', READY, '4k'),
    ]),
    (5, '药物和最终靶点', [
        ('a', '5a_antigen_funnel', '抗原级完整漏斗阶梯：213 个抗原各占一列、按淘汰门排序，10 行门（全部、C1–C5、MR、coloc、单一修复方向、方向匹配药物），格色 = 通过或在此淘汰，最终靶点着色，右列剩余数（148 → 9 → 7 → 5 → 2）', '1x4', 1, 1,
         'R03/target_funnel_final_20261119/antigen_gate_matrix.csv、funnel_counts.csv、near_misses.csv', '#26 waffle()；#20 条码（HH/Fig4/G.R）', READY, '5a'),
        ('b', '5b_program_flow', '程序级流向三轴冲积：程序（6,849）→ 所属抗原 → 最后通过的门', '1x2', 2, 1,
         'R03/target_funnel_final_20261119/antigen_gate_matrix.csv；OPT/all_program_evidence.csv', '#27 AD v3 图 5c（ad_omics_v3/delivery/figures/source/render_current.R）；dribbon()', READY, '5b'),
        ('c', '5c_evidence_wheel', '主视觉。证据轮：扇区 = 两个最终靶点，角向列 = 程序（CD71 4 个、CD119 6 个），8 环 = C1–C5、MR、coloc、方向匹配药物；外柱 = 满足的门数，全满程序突出；最外细刻度标 BIG40、反向 MR、平台一致（注释，不是门）', '2x2', 2, 3,
         'R03/target_funnel_final_20261119/final_targets.csv、final_target_mr_coloc_pairs.csv、final_target_drug_matches.csv', '#28 PeriScope 证据轮（evidence_wheel.py、rings_grid.R）；dring()', READY, '5c'),
        ('d', '5d_drug_matrix', '方向匹配药物属性符号矩阵（6 个：EMAPALUMAB、AMG-811、FONTOLIZUMAB、A27.15、E2.3、PPMX-T003）：最高期别（位置与颜色）、作用类型（形状）、途径（直接、配体途径）、与模型修复方向一致（实心），按靶点分组', '1x1', 3, 1,
         'R03/target_funnel_final_20261119/final_target_drug_matches.csv；R03/drug_evidence_expanded_20261119/drug_records.csv.gz', '#29 CREED paper2 1b', READY, '5d'),
        ('e', '5e_drug_ledger', '药物记录账本（原 e、f 合并）：两个靶点的 56 条药物记录一格一条，行 = 途径（直接、配体途径、复合物伙伴、基因沉默），列组 = 相对模型方向（匹配、相反、细胞杀伤或方向未知），填色 = ChEMBL 模态', '1x1', 3, 2,
         'R03/drug_evidence_expanded_20261119/drug_records.csv.gz、antigen_drug_summary.csv', '#26 waffle()；hh S7a 单位堆叠（HH/S7/A.R）', READY, '5e + 5f'),
        ('f', '5f_gene_counts', '程序基因数热表：行 = 靶点、定义、方向，列 = 8 条路线，格内印基因数，行与列合计点', '1x2', 4, 1,
         'OPT/all_program_evidence.csv（n_genes、genes）', 'hh S9a 热表（HH/S9/A.R）', READY, '5g'),
        ('g', '5g_gene_upset', '单位 UpSet：10 个最终程序的基因并集，一格一个基因，每 10 格一列；集合大小也用方块', '2x2', 4, 3,
         'OPT/all_program_evidence.csv（genes）', '#30 hh 图 5f（HH/Fig5/F.R）', READY, '5h'),
        ('h', '5h_shared_genes', '共享基因清单：共享基因 × 程序成员瓦片（行 = 基因，列 = 程序，按靶点分块）', '1x2', 5, 1,
         'OPT/all_program_evidence.csv（genes）', '#20 条码热图', READY, '5i'),
    ]),
]


def grid_map(panels):
    cells = [['.'] * 4 for _ in range(5)]
    for p in panels:
        r, c = (int(v) for v in p['block'].split('x'))
        for i in range(p['row'], p['row'] + r):
            for j in range(p['col'], p['col'] + c):
                assert cells[i - 1][j - 1] == '.', (p['id'], i, j)
                cells[i - 1][j - 1] = p['letter'] if p['status'] == READY else p['letter'].upper() + '*'
    return cells


def units(b):
    r, c = (int(v) for v in b.split('x'))
    return r * c


md = ['# PeriScope 新主图版面（figures_v2_20261119）', '',
      '2026-11-19。本文件和 `layout_fig1.json` 到 `layout_fig5.json` 由 `gen_layouts.py` 一并生成，改版面时改那个脚本再重新生成。',
      '', '## 网格和规则', '',
      '- 单位格 u = 41.25 mm，间隔 g = 3 mm。页宽 174 mm = 4u + 3g；最多 5 行，页高 218.25 mm。',
      '- "块"写作 行×列（高×宽，单位格）：1×1 为 41.25 × 41.25 mm，1×2 为 85.5 × 41.25（宽 × 高），2×1 为 41.25 × 85.5，1×3 为 129.75 × 41.25，1×4 为 174 × 41.25，2×2 为 85.5 × 85.5。允许的块：1×1、1×2、1×3、1×4、2×2、2×3、2×4 及其转置。',
      '- "位置"是块左上角所在的格（行，列），从 1 数起。每张图 20 个单位格全部用满。',
      '- 面积按信息密度分：多编码、多单位的组合图给 2 到 4 格，单一简单图给 1 格。',
      '- 子图左上角 3.5 × 3.2 mm 留空，拼页时在这里印 8 pt 粗体小写字母。',
      '- 待补内容（消融）只画细框空块，不印任何文字，也不印字母；后面的字母照常保留。',
      '- 字体、字号、线宽、配色、检查和拼页工具见 `_system/README.md`。',
      '', '缩写：' + '；'.join(f'{k} = `{v}`' for k, v in ABBR.items()) + '。', '']

for n, title, rows in FIGS:
    panels = [dict(letter=l, id=i, content=c, block=b, row=r, col=co, source=s, reuse=ru, status=st, plan=pl,
                   file=f'fig{n}/panels/{i}_{b}.pdf') for l, i, c, b, r, co, s, ru, st, pl in rows]
    cells = grid_map(panels)
    assert all(x != '.' for row in cells for x in row), f'figure {n} leaves empty cells'
    assert sum(units(p['block']) for p in panels) == 20
    def expand(s):
        for k, v in ABBR.items():
            s = s.replace(k + '/', v.split(' (')[0] + '/')
        return s
    js = {'figure': n, 'title': title, 'grid': {'u_mm': U, 'g_mm': G, 'cols': 4, 'rows': 5}, 'letter_on_blank': False,
          'panels': [dict(p, source_paths=expand(p['source'])) for p in panels]}
    with open(os.path.join(OUT, f'layout_fig{n}.json'), 'w') as f:
        json.dump(js, f, ensure_ascii=False, indent=1)
    big = max(panels, key=lambda p: units(p['block']))
    md += [f'## 图 {n}　{title}（{len(panels)} 个子图）', '', '```text', '     列1  列2  列3  列4']
    md += [f'行{i + 1}  ' + ''.join(f'{x:<5}' for x in row).rstrip() for i, row in enumerate(cells)]
    md += ['```', '']
    if any(p['status'] == BLANK for p in panels):
        md += ['大写加星号的格是待补空块。', '']
    md += ['| 字母 | 内容 | 块（行×列） | 位置（行，列） | 数据来源 | 借用设计 | 状态 | 原方案 |', '|---|---|---|---|---|---|---|---|']
    for p in panels:
        md.append(f'| {p["letter"]} | {p["content"]} | {p["block"].replace("x", "×")} | {p["row"]}，{p["col"]} | {p["source"]} | {p["reuse"]} | {p["status"]} | {p["plan"]} |')
    md.append('')

md += ['## 相对 PLAN.md 的改动', '',
       '内容全部保留，只做合并、拆分和重排，以便按信息密度分面积、每张图正好用满 20 格。',
       '',
       '- 图 1：原 g（留出 CITE-seq 评估示意）并入 a 的右端评估支路；原 1E 标记点热图并入 f（冲积）右侧；字母按阅读顺序重排（原 c→e、d→d、e→f、f→c、h→g、i→h、k→i、j→j、l→k）。',
       '- 图 2：原 h（种子稳定性）提前为 d，原 d（外部队列）为 e，原 i（读出优化）为 f；三个消融（原 e、f、g）放在第 4 行，暂为空块 g、h、i。b 瓦片内的逐种子小柱改为小点，因为比较方法只有 3 个种子，小柱会成为"3 根柱"的柱图。',
       '- 图 3：原 b（基因响应点热图）和 c（剂量—响应卡）合成一张五卡组合图 b；原 f、g 合成 f（种子热条加匹配零分布），原 h、i 合成 g（留一种子方块加三队列供体效应），f 与 g 同处第 3 行、共用行序；原 d 为 c，e 不变，原 j、k、l、m、n 依次为 d、i、j、h、k。行集合为两个最终靶点的 10 个程序。',
       '- 图 4：SuSiE-coloc 已完成，g 不再占位；可选的区域关联图 k 以 1×1 放入，用最终靶点之一。',
       '- 图 5：原 e（药物筛除账本）和 f（药物模态）合成一个账本 e（一格一条记录，填色为模态），原 g、h、i 依次为 f、g、h。按用户要求不放"药物方向与模型方向示意"和"最终靶点总结卡"。',
       '- 全部图：不使用 2026-11-19 回滚的种子 52–61；完整模型种子为 42–51，模型方向取确认种子 47–51；不画任何无模型基线（donor 均值、pseudobulk PD 减对照）。', '']

with open(os.path.join(OUT, 'LAYOUT.md'), 'w') as f:
    f.write('\n'.join(md) + '\n')
print('wrote layout_fig1..5.json and LAYOUT.md to', OUT)

from pathlib import Path
import pandas as pd,numpy as np,json,math,hashlib,shutil
W=Path(__file__).resolve().parent;O=W.parents[2]/'delivery/blood_brain_MR_20261003'
d=pd.read_csv(O/'all_brain_MR_estimates.csv');v=d[d.eligible].copy();b=v['sign']*v.beta_outcome/v.beta_exposure;se=v.se_outcome/abs(v.beta_exposure);p=np.array([math.erfc(abs(x)/math.sqrt(2)) for x in b/se])
assert np.allclose(b,v.beta_MR,atol=1e-12) and np.allclose(se,v.se_MR,atol=1e-12) and np.allclose(p,v.p_MR,atol=1e-12)
maxq=0.
for name,g in v.groupby('analysis_set'):
 a=g.p_MR.to_numpy();order=np.argsort(a);q=np.minimum.accumulate((a[order]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1].clip(0,1);maxq=max(maxq,float(np.max(abs(q-g.q_BH.to_numpy()[order]))))
assert maxq<1e-12
m=json.loads((O/'brain_data_acquisition.json').read_text());assert m['complete'] and m['n_AGES_pairs']==224 and m['n_outcomes']==16
primary=v[v.analysis_set.eq('independent_AGES_primary')];assert len(primary)==192 and primary.gene.nunique()==12
check=dict(status='passed',n_acquired_pairs=len(d),eligible_estimates=len(v),primary_independent_tests=len(primary),unique_primary_genes=12,brain_outcomes=16,max_beta_error=float(np.max(abs(b-v.beta_MR))),max_se_error=float(np.max(abs(se-v.se_MR))),max_p_error=float(np.max(abs(p-v.p_MR))),max_BH_error=maxq,checks=['brain/exposure position and rsid explicitly checked in same genome build','all sixteen outcomes retained regardless of P','original exposure beta and EA matched to published supplement','all HTTP ranges and BGZF blocks validated','MR independently recomputed outside R','independent serum and overlapping UKB plasma analyses separated'])
(O/'verification.json').write_text(json.dumps(check,indent=2))
uk=v[v.analysis_set.eq('UKB_overlap_sensitivity')];hits=uk[uk.q_BH.lt(.05)]
lines=['血端蛋白 → 脑表型 MR：数据和首轮结果','',
'脑表型固定为双侧尾状核、壳核、苍白球、丘脑的体积与T2*，共16项。GWAS Catalog全基因组发现队列文件：体积N=22128，T2* N=20043；没有用网页的显著关联名单替代完整统计量。',
'成功取得865组蛋白测定×脑表型工具变异统计量：AGES224组、UKB641组。所有16个结局查询完整结束，无网络错误被记成阴性或缺失。少数UKB变异没有精确匹配，见提取审计。',
'独立血端：AGES冰岛血清，14项测定覆盖12个当前候选；每个候选按血端F统计量选一项主测定，12×16=192项独立两样本MR。另两项测定保留为敏感性，回文等位基因缺少脑端频率时不强行定向。',
'UKB血浆：43项既有cis工具测定；通过等位基因匹配的560项MR单独分析。血端与MRI样本可能重叠，不计为独立复现。',
'deCODE冰岛血浆：已下载公开补充表，23项测定覆盖22个候选。补充表beta经过四舍五入且缺少完整SE；现阶段只列为可扩展数据源，不混入完整精度MR。官方完整数据入口需要填写姓名、邮箱、机构及验证码，尚未提交表单。','',
'首轮结果',
f'AGES独立主分析：192项，名义P<0.05共{int((primary.p_MR<.05).sum())}项；FDR<0.05共{int((primary.q_BH<.05).sum())}项。',
f'UKB分析：560项，名义P<0.05共{int((uk.p_MR<.05).sum())}项；FDR<0.05共{len(hits)}项，涉及TFRC/CD71、CD40、CD22。']
for _,r in hits.iterrows():lines.append(f"{r.gene} | {r.outcome_trait} | beta={r.beta_MR:.5g}, SE={r.se_MR:.5g}, P={r.p_MR:.4g}, FDR={r.q_BH:.4g}")
lines+=['','CD22对应的是右侧苍白球体积关联。它提供了可继续检验的血端蛋白—脑表型遗传线索；当前仍需区域共定位和独立数据支持，不能由本轮单SNP MR单独确定外周介导路径。T2*或体积变化的方向不直接等同于治疗获益。','',
'统计方法：单cis工具Wald比值、双侧检验、一阶标准误；另存二阶delta-method标准误。主分析、替代测定和UKB重叠样本分析分别做BH校正。AGES效应单位为原文Box-Cox转换后的蛋白分析单位，脑表型单位为标准化表型SD。不是PD患病风险OR。','',
'复现：来源、坐标、等位基因、beta/SE和逐项匹配记录见CSV与JSON。计算脚本和协议在reproduce/；原始公开表及逐块下载缓存位于过程目录blood_brain_MR_20261003。无需重新下载即可从本目录的combined_exposure_instruments.csv与brain_sentinel_statistics.csv运行R脚本复现首轮MR。所有结局和阴性估计均保留。']
(O/'MR数据与结果说明.txt').write_text('\n'.join(lines)+'\n')
rep=O/'reproduce';rep.mkdir(exist_ok=True)
for n in ['run_brain_mr.R','analysis_protocol.txt']:
 shutil.copy2(W/n,rep/n)
(O/'manifest.json').write_text(json.dumps({str(p.relative_to(O)):hashlib.sha256(p.read_bytes()).hexdigest() for p in O.rglob('*') if p.is_file() and p.name!='manifest.json'},indent=2))
print(check)

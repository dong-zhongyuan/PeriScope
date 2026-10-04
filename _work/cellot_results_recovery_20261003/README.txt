CellOT 结果复现与恢复

训练和评价函数沿用误删前版本；其他方法与 MR 不参与此流程。
标准入口：run_verified_recipe.py。服务器 Python：/public/home/mengxl/dzy/envs/pd_bbm/bin/python。
运行目录：/public/home/mengxl/dzy/pd_product/_work/cellot_results_recovery_20261003。
先运行 run_verified_recipe.py，旧值校验通过后运行 install_recovery.py；本地运行 complete_and_sync.py 可等待校验、安装并核对下载文件哈希。
迁移机器时，统一修改脚本中的服务器项目、资产与本地项目根路径，并使用 recovery_environment.json 指定的依赖版本及相同输入文件。

已确认的原数值设置：
原六条轴的 PCA：BLAS 双线程；pDC 两条新增轴的 PCA：BLAS 单线程。均为 PCA64、randomized、seed42。
原六条轴的 seed42 训练：CPU，Torch 双线程，BLAS 双线程。
其余18个种子×轴：CPU，Torch 单线程，BLAS 单线程。
残差读出：BLAS 单线程；原预测转 float64 后加原训练重构残差。
猪预测：BLAS 双线程，Torch 双线程，原 CUDA 推理路径。人类训练权重不使用猪脑结果调节。
核对标准见 expected_results 中找回的原始数值以及 validate_recovery.py；不按测试表现挑选新参数。

最终文件：
benchmark_20261002/models、results、readout_repair/results：权重、原预测、残差读出。
benchmark_20261002/cellot_projection：已验证 PCA，可直接保留复用。
benchmark_20261002/cellot_restored：旧数值、逐项校验、运行设置及输入哈希。
pig_external_validation_20261003/cellot_restored：CellOT 猪预测、置换分布和各项明细。

gpu_seed42.py、finish_recovery.py 属于恢复期间的运行设置排查工具，不是标准复现入口。

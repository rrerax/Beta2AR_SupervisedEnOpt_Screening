# Preprint Outline — Supervised EnOpt Virtual Screening (β2-AR), with an FtsZ Transfer Case Study

> Status: v0 draft outline · 2026-09-05
> 用法:本文件是把两个仓库(Beta2AR_EnOpt_GPCR_Screening 主线 + FtsZ_Coumarin_HTVS_Reconstruction case study)写成一篇方法学 pre-print 的"施工图"。每一节都标注了:写什么、从哪个仓库文件取素材、哪些素材还需要新做。

## 0. 决策日志 (Decisions)

- [x] 组织方式:**单篇 preprint(已确认 2026-09-05)**。主线 = β2-AR 监督 EnOpt;FtsZ 作为"同一框架在无实验标签场景的可迁移性"case study 小节。决策背景:作者硕士答辩已完成,FtsZ 无需独立可引用成果;两个代码仓库保持独立并互链,论文层面合一。
- [ ] 目标平台:建议先投 **bioRxiv**(bioinformatics/methods 板块)抢占时间戳,再视反馈投期刊(JCIM / J. Chem. Inf. Model. 这类纯计算方法学期刊最对口)。
- [ ] 托管位置:正式手稿建议新建独立仓库或放 Beta2AR 仓库 `paper/` 目录;本文档先放在本仓库 `docs/` 供两仓库共用引用。
- [ ] 语言:先英文成稿(目标期刊语言),中文版留给答辩/组会材料。

## 1. 论文定位

### 1.1 一句话主线
**把"多构象对接 + 用已知活性学习分数组合"(监督 EnOpt)做成一整套可复现工作流,并诚实回答三个问题:它比简单平均好在哪、系统性偏差(大分子虚高)怎么修、采样质量够不够。**

### 1.2 支持的主张 (claims,每一条都要有数据支撑)
1. 对 β2-AR,用 5 个活性/非活性构象做 ensemble docking 能稳定重现(Stage 1)。
2. 监督 EnOpt 的价值集中在**榜单头部**(top-1% EF ~2× oriented mean),而不是整体 AUROC;在 matched-decoy 硬测试上与正确朝向的简单平均统计持平(0.696 vs 0.719,CI 重叠)。
3. 追加 per-heavy-atom 分数可以**同时**修掉"大分子虚高"偏差并提高判别力(0.696→0.751,EF1 3.9→7.3)。
4. 提高 exhaustiveness(ex4→ex8)不改变结论(CI 重叠):瓶颈在受体/姿势建模,不在采样。
5. 无标签小规模场景(旧数据重构)仍可从 ensemble 重排获益(排序稳定性),但**不能宣称活性**。

### 1.3 明确不主张 (non-claims)
- 不声称任何分子有实验活性;榜单只是采购/验证候选。
- 不声称监督 EnOpt 是普适优于 consensus 的方法(证据:硬测试两者统计持平)。
- FtsZ 部分不声称是完整 EnOpt 模型(无标签,仅为 EnOpt-style consensus)。

## 2. 备选标题

- EN1: *Learning to rank docking-score ensembles: a reproducible supervised EnOpt screen of the β2-adrenergic receptor and a label-free transfer case on FtsZ*
- EN2: *What machine learning does and does not add to ensemble docking: a β2-AR case study with DUD-E validation*
- EN3: *Reproducible ensemble-docking screening with supervised reranking: honest benchmarks for β2-AR and a legacy FtsZ screen*
- CN(参考): *多构象对接的监督式重排:β2-肾上腺素受体的可复现筛选与 FtsZ 可迁移案例*

## 3. 手稿章节骨架

### Abstract (DRAFT v0,约 200 词,先改后定稿)

High-throughput virtual screening against a single rigid receptor structure is cheap but unreliable: one scoring run cannot capture the conformational plasticity of G-protein-coupled receptors, and raw docking scores carry systematic biases such as molecular-size dependence. The ensemble-optimizer (EnOpt) concept addresses the first problem by combining per-conformation docking scores into one ranking, but its value depends on how scores are combined. We present a reproducible two-stage workflow that applies supervised EnOpt reranking to the β2-adrenergic receptor (β2-AR). 29,865 library molecules were docked against five active/inactive β2-AR conformations (149,325 docking runs, AutoDock Vina 1.2). An XGBoost model trained on known β2-AR actives (48, later expanded to 206 literature actives) and validated against 2,978 property-matched DUD-E decoys separated actives from decoys with out-of-fold AUROC 0.696 (95% CI 0.655-0.737), statistically indistinguishable from a correctly oriented five-pocket average but with roughly twofold top-1% enrichment. Appending per-heavy-atom-normalized scores raised AUROC to 0.751 (0.712-0.791) and removed the large-molecule ranking bias. Re-docking the validation set at exhaustiveness 8 left all conclusions unchanged (confidence intervals overlap), showing the remaining bottleneck is receptor and pose modelling rather than docking sampling. A chemical review of the top-200 shortlist shows the recommended candidates are diverse and novel (170 scaffolds; 192/200 have max-Tanimoto < 0.4 to training actives). As a transfer case study, the same ensemble framework was applied to a legacy coumarin/FtsZ screen reconstructed from thesis spreadsheets without experimental labels, improving ranking stability but not permitting activity claims. All code, score tables, trained models, and reports are provided.

### 1 Introduction (素材: 两个 README + docs/method_notes.md)
- 背景:VS 常用单构象刚性受体;GPCR 的 TM6 激活运动改变口袋 → 需要构象集合(引 EnOpt 论文 Bhatt et al. 2024;顺带 β2-AR 有丰富晶体结构)。
- 问题1:ensemble docking 产生"每个分子一组分数",怎么合成?(简单平均、best score、训练模型)。
- 问题2:评分函数系统偏差——Vina 分数与分子大小相关,大分子虚高。
- 问题3:ML 在 VS 里怎么用才可信——OOF 防自评、matched decoys、CI、失败案例报告。
- 目标:给出可复现 + 诚实评估的模板,并在第二个靶点(无标签旧数据)演示框架可迁移。

### 2 Methods (素材: scripts/ 01-11 + revalidation_ex8/ + FtsZ scripts/)
- 2.1 通用框架定义(建议画成 Fig 1):
  `输入库 → 配体/受体准备 → 多构象对接 → 分子×构象分数矩阵 → 合成器(consensus 或 supervised) → 验证 → 候选短名单`
- 2.2 Case A(β2-AR):
  - 受体:2RH1/5D5A(inactive)+ 4LDE/3SN6/4LDL(active),5 构象。
  - 库:ChEMBL37 子集 30,000;准备成功 29,865;Vina 1.2 ex4。
  - 正例:48 个库内已知活性 → 扩充到 206(ChEMBL pChEMBL≥6、置信≥8)。
  - 负例:DUD-E decoys 2,978(性质匹配)。
  - 特征/模型:XGBoost 3 折 OOF;raw 5 特征;raw+per-ha 10 特征;raw+ha+HAC 11。
  - 稳健性:ex8(exhaustiveness 8, num_modes 20)重对接验证集与 top-200。
- 2.3 Case B(FtsZ):
  - 旧工作簿重构(列含义按内容推断)、RDKit 去盐/规范化、双口袋(ezPocket 中心/残基锚)、Vina rerun、ANM 5 构象 ensemble、EnOpt-style 权重(无标签→只做秩一致性共识)。
- 2.4 验证协议:3 折 OOF、AUROC+95%CI、EF1/5、化学评审(Butina/Bemis-Murcko/Tanimoto 新颖性/PAINS-Brenk)。

### 3 Results (素材: README Stage 1-5 表 + results/ 各目录)
- 3.1 β2-AR 数据集与对接(Stage1 数字)。
- 3.2 Stage2:48 正例时代 OOF 0.6647(说明当时 negative 假设的问题)→ 被 Stage3 取代。
- 3.3 Stage3 DUD-E 硬测试:0.696 vs oriented mean 0.719(CI 重叠);合并账目 0.718 vs 0.671;top-1% EF 3.88 vs 1.94。
- 3.4 Stage4 per-heavy-atom:0.696→0.751(EF1 7.28);尺寸偏差修复(MW top200 393→358)。
- 3.5 Stage4a 回顾与 top-200 化学评审:库里仅 50 个已知活性 → 回顾无统计功效;评审给出 170 scaffolds/192 novel/75 PAINS。
- 3.6 Stage5 ex8:验证集 0.654(0.612-0.695)/0.733(0.692-0.772),与 ex4 CI 重叠;top-200 排序稳定(Spearman 0.87-0.91,89/100 重合);CHEMBL776 仍垫底 → 系统性欠分案例。
- 3.7 FtsZ case:90 结构 450 dock;legacy vs ensemble 相关(BP2 Spearman ~0.80;BP1 ~0.53);top weighted hits 展示;无标签边界。

### 4 Discussion & Limitations (素材: 各 docs/*_notes.md 的 Honest reading)
- 4.1 何时 ML 有用:榜单头部富集、整库账目;何时没用:matched-decoy 整体 AUROC。
- 4.2 per-ha 是便宜且可复用的"偏差修复"特征;教训:特征要有化学意义。
- 4.3 采样 vs 模型:ex8 结果指向 receptor flexibility/IFP/pose 校验是下一瓶颈。
- 4.4 失败案例透明化(CHEMBL776)是加分项,建议保留为 figure/SI。
- 4.5 局限:无湿实验;decoy 匹配依赖性质;ChEMBL 回顾无统计功效;FtsZ 无标签 + ANM 构象 + 重建 box。

### 5 Data & Code Availability
- 归属与溯源:FtsZ 重构部分是作者硕士项目成果、早于本工作;在 Acknowledgement / Methods 中写明 lineage("the ensemble reranking idea was piloted in the author's FtsZ reconstruction, refs/URL")。β2-AR 为本工作主线。
- Beta2AR:数据全部 ChEMBL-derived 可公开;30k 库、分数表、模型、榜单在 `results/`、`data/`(models 是 XGBoost pickle)。
- FtsZ:重构表与结果公开;原始 thesis 工作簿与受体 PDB 在本地 `raw_inputs_staging/`(不随仓库发布),数据可用性说明见 FtsZ `docs/data_availability.md`。
- 环境:`environment-vina.yml` / `requirements.txt`。

## 4. 图/表计划(映射到仓库文件)

| 编号 | 内容 | 素材来源(现有) | 需要新做 |
|---|---|---|---|
| Fig 1 | 通用框架示意(两靶点共用) | 无 | 新画(建议矢量图;风格参考 FtsZ workflow_summary.png) |
| Fig 2 | β2-AR 工作流/数据规模总览 | FtsZ `results/figures/workflow_summary.png`(风格参考) | 新做一张 β2-AR 版 |
| Fig 3 | 构象分数分布与 top hits | Beta2AR `results/figures/score_distributions.png`、`enopt_weighted_top_hits.png`;`supervised_enopt/fig1_score_comparison.png`/fig2/fig3 | 按投稿尺寸重排 |
| Fig 4 | DUD-E 验证结果(ROC/EF/CI) | `supervised_enopt_decoy/fig_decoy_validation.png` | 可选合并成投稿版 |
| Fig 5 | per-heavy-atom 特征对比 | `feature_experiment_heavy_atom/fig_feature_comparison.png` | 高亮 0.696→0.751 |
| Fig 6 | ex4 vs ex8 稳健性 | `ex8_rerank/fig_ex4_vs_ex8.png` | 补充 CHEMBL776 标注 |
| Fig 7 | top-200 化学评审(MW/HAC、新颖性、PAINS) | `review_top200_ha/fig_top200_properties_vs_actives.png` | 可加化学空间 t-SNE |
| Fig 8 | FtsZ case:rerun 一致性 + ensemble 重排 | FtsZ `results/figures/vina_rerun_vs_original_scores.png`、`enopt_style_reranking_vs_legacy.png`、`ensemble_score_matrix_heatmap.png` | 裁剪合并 |
| Table 1 | 两靶点数据集总览 | 两 README Summary 表 | 汇总成一张 |
| Table 2 | β2-AR 验证主结果(raw / raw+ha / ex4 / ex8 + EF) | README Stage 3-5 表 | 汇总 |
| Table 3 | FtsZ legacy-vs-rerun / ensemble 相关 | FtsZ `results/tables/enopt_style_correlations.csv`、`vina_rerun_correlations.csv` | 汇总 |
| Table 4(SI) | 完整榜单/评审/模型卡索引 | 见 §6 SI 清单 | 无 |

## 5. 数字速查表(全篇必须一致,写稿时以此为准)

| 量 | β2-AR 值 | 备注 |
|---|---|---|
| 库规模 | 30,000 → 29,865 可对接 | 149,325 rows |
| 构象 | 5(2 inactive + 3 active) | 2RH1/5D5A/4LDE/3SN6/4LDL |
| Stage2 正例 | 48(库内) | OOF 0.6647(基线未校正,已被取代) |
| Stage3 正/负 | 206 / 2,978 DUD-E decoys | 硬测试 |
| Stage3 OOF AUROC | 0.696(0.655-0.737) | EF1 3.88 / EF5 2.72 |
| oriented mean(硬测试) | 0.719(0.679-0.760) | 与 EnOpt CI 重叠 |
| 合并账目 AUROC | EnOpt 0.718 / mean 0.671 | 全库+decoys+新活性 |
| raw+per-ha AUROC | 0.751(0.712-0.791) | EF1 7.28 / EF5 4.76 |
| top200 MW(Stage3→raw+ha) | 393 → 358 | 活性窗口中位 ~362 |
| ex8 验证(raw/raw+ha) | 0.654(0.612-0.695)/0.733(0.692-0.772) | 与 ex4 CI 重叠 |
| ex4 top100 留在 ex8 top100 | 89/100(best) | Spearman 0.87-0.91 |
| CHEMBL776 | rank 27,659/29,865;ex8 201/201 | 失败案例 |

| 量 | FtsZ 值 | 备注 |
|---|---|---|
| 记录/结构 | 204 method records → 90 unique | BP1 55 / BP2 72 化合物 |
| ensemble | 5 构象(ANM) | 450/450 docks |
| legacy vs ensemble Spearman | BP1 ~0.53;BP2 ~0.80 | enopt_style_correlations.csv |
| 标签 | 无 | EnOpt-style,非 supervised |

## 6. 参考文献(起步集,写作时用 citation-management skill 核验)

1. Bhatt NM, Wang A, Durrant JD. *Sci Rep* 2024;14:20722 (EnOpt 概念)。
2. Mysinger MM, et al. *J Med Chem* 2012;55:6582-6594 (DUD-E)。
3. Trott O, Olson AJ. *J Comput Chem* 2010;31:455-461 (AutoDock Vina)。
4. Eberhardt J, et al. *J Chem Inf Model* 2021;61:3891-3898 (Vina 1.2)。
5. Mendez D, et al. *Nucleic Acids Res* 2019;47:D930-D940 (ChEMBL)。
6. Bemis GW, Murcko MA. *J Med Chem* 1996;39:2887-2893 (scaffolds)。
7. Baell JB, Holloway GA. *J Med Chem* 2010;53:2719-2740 (PAINS)。
8. RDKit: Landrum G, rdkit.org(ECFP4/Butina/去盐)按需补 Rogers & Hahn 2010。

## 7. 可选:两篇拆分方案(备选,若审稿人/导师希望更聚焦)

- **Paper A(β2-AR 主线)**:上述 3.1-3.6 + Discussion,删掉 FtsZ 小节 → 更"深"。
- **Paper B(工具/流程)**:把通用框架 + 两案例写成"可复现 ensemble-screening 模板 + 诚实评估清单",更像 software/method note。
- 现状决策:答辩已完成 → **已锁定单篇**,本备选仅在 FtsZ 未来想单独投 reproducibility note 时启用。
- 决策标准(备忘):若 β2-AR 后续补了 IFP/柔性侧链实验并冲到 AUROC 0.8,可再考虑拆 Paper A。

## 8. 下一步 checklist

- [ ] 定标题 + 期刊 + 作者分工(需你拍板)
- [ ] 用 Abstract DRAFT v0 起稿;Intro 起草可用本地 research-writing / nature-writing skill
- [ ] 新做 Fig 1(框架示意)+ Fig 2(β2-AR 总览),其余图先复用
- [ ] 汇总 Table 1/2/3(数据已在速查表)
- [ ] 把 outline 拆成 `paper/manuscript.md` 逐节填充
- [ ] 数据可用性:确认 30k 库来源文件可公开、模型 pkl 版本说明

## 9. 需要你拍板的决定

1. ~~单篇还是两篇?~~ **已定:单篇(2026-09-05,硕士答辩已完成)**。
2. 标题选 EN1/EN2/EN3 还是新拟?【待定】
3. 手稿仓库放哪:Beta2AR 仓库 `paper/`、FtsZ 仓库、还是新建独立仓库?
4. 第一作者与贡献者清单(决定"我们"的表述和致谢)。
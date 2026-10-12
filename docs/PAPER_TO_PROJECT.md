# 原论文与本项目的对应关系

依据用户提供的 25 页 PDF *Rethinking Molecular Text Representations for LLMs: An Empirical Study*（arXiv:2606.03057v1）。下表页码是 PDF 页码，引文为原文实际文字。

**关于作者代码：本 PDF 全文做 URL 扫描后只出现一个链接 `https://huggingface.co/mistralai/`（模型出处），未发现作者的 benchmark 源码仓库链接。** 因此这里映射的是论文文字与本项目实际源码，不声称检查过作者的运行源码。

原文规模（摘要与 §1 实际文字）："nine representations and eight chemical tasks"、"We benchmark 16 LLMs across five model families"、"the 16 models with the 9 molecular representations across the 8 tasks"。

| 原文位置 | 论文方法 | 本项目代码与差异 | 采用/修改/未实现 |
| --- | --- | --- | --- |
| 摘要 / §1 | 9 种表示、8 项任务、16 个模型 | `bmb/reps.py`、`reps_ext.py`：8 种输入（SMILES、SELFIES、Our Code、img、Graph2D、Graph3D、CML、MolJSON）；5 项任务；1 个模型 | 修改：未涵盖原文全部 SMILES 变体、IUPAC、InChI；只测一个服务模型别名 |
| p15 Appendix C | "We use the test split as the source of molecules"；从 ~3.3K 取 "a stratified subset of 200 molecules" | `dataset.py/read_source_rows` 读 ChEBI-20 **再发布**的 validation + test 再合并，由项目重划 dev/test；候选池 431 | 修改：不是官方 test 划分 |
| p15 Appendix C | "molecular weight ≥300 Da and at least 2 rings" | `dataset.py/apply_complexity_filter`：沿用该门槛；具体用 `Descriptors.MolWt` 与 `CalcNumRings` 是**本项目口径选择**，论文未指定实现 | 采用门槛 + 项目自定义实现 |
| p15 Appendix C | 另用 ZINC250K 构造 tautomer / protonation 配对 | 本项目不使用 ZINC250K，也没有这两项任务 | 未实现 |
| p15 Appendix D | Atom Counting：逐次问 C/N/O/S/F/Cl 中"a specified element"，exact match | `labels.py/atom_count` 仍是 RDKit 原子迭代；`tasks.py/AtomCountTask.messages` **一次问四个元素**，四值全对才得分 | 修改：问法与元素集合（用 S 而非 F，理由见 `labels.py` 注释中的分布统计） |
| p15 Appendix D | Functional Group：5 类（aldehyde, ester, halide, primary amine, sulfonamide），逐次 Yes/No，macro-F1 | `labels.py/functional_group_labels` 用 SMARTS 匹配；一次问 **13 类**，主分数固定 **10 类**；13 类中 aldehyde / nitrile / nitro 被问到但不进主分 | 修改：类别集合与问法 |
| p15 Appendix D | Molecular Property Estimation：LogP、TPSA、HBD、HBA，"evaluated by exact match accuracy" | `DescriptorTask` 本轮只问 HBD；`QedTask` 是本项目**新增**的 QED 估计，用 MAE | 部分采用 + 新增：不能说论文采用了 QED |
| p15 Appendix D | Molecule Retrieval：四选一，top-1 | 无 | 未实现 |
| p15 Appendix D | Isomer Discrimination：500 pairs，含同分子异遍历、**立体异构**、替换三类 | `tasks.py/IdentityTask/build_pairs`：200 对，无立体化学（池本身无立体信息），正例重排/重遍历、负例取不同结构 | 修改：不是原任务复现 |
| p15–16 Appendix D | Tautomer Recognition（250 对）、Protonation State Recognition（210 对） | 无 | 未实现 |
| p15–16 Appendix D | Caption-to-Molecule Generation：validity / exact match / Tanimoto / FCD | 无 | 未实现 |
| p16–17 Appendix E | 统一任务 prompt，"reason step by step, and put your final answer within \boxed{}" | `tasks.py/system_message` 要求 "a single compact JSON object and nothing else"；无逐步推理要求 | 修改：输出格式与推理要求都不同 |
| p3 §3 | "logP and TPSA are regression tasks using Spearman correlation ρ, while hydrogen-bond donor/acceptor counts use exact-match accuracy" | 本项目 QED 明确用 MAE、HBD 用精确匹配，并在报告中写明自己的口径 | **原文内部不一致**：§3 对 LogP/TPSA 说 Spearman，Appendix D 却说 property estimation 整体 "evaluated by exact match accuracy"。本项目不继承任一说法，只报告自己实际使用的指标 |
| p4 §4 | "a paired bootstrap test (Koehn, 2004) with 10,000 resamples between the per-instance scores of every other representation and the top scorer" | 项目使用冻结配对、重采样与 Holm 校正，具体以 `scoring.py`、`build_analysis.py` 的当前版本为准 | 部分采用：**不能声称完全照搬**原文统计方法；Identity 的确认性检验未获验证 |
| p7 §5.4 | 用 Gemini 3 Flash 作为 judge 评判生成错误模式 | 无生成任务、无 LLM judge | 未实现 |
| p7–9 §6，附录 H/I | tokenization 审计、hidden-state 线性 probe、attention 分析 | 只记录 API 用量，没有隐藏状态、没有训练 probe | 未实现：用量记录不等同机制研究 |
| p9 Limitations | 只比较文本，提出多模态为未来工作 | 图片与 2D/3D 文本是本项目扩展；但 MOL 3D 文本不等于完整三维推理能力证据 | 扩展 |

## 需要向审核者说明的三点

1. **Identity 不等于原论文的 Isomer Discrimination。** 原文 500 对含立体异构样本，本项目池无立体化学，且历史构造中三表示与五表示的视图规则不同（见 `docs/AUDIT.md`），所以 Identity 成绩只作描述性比较。
2. **数据划分不是原划分。** 本项目合并 ChEBI-20 再发布的 validation 与 test 后自己重划，dev 用于调方案、test 用于报告；不能直接说"用了论文的 200 个样本"。
3. **分数不可与原文横向比较。** 数据、标签、问法、表示实现、模型都不同；我们的 CML/MolJSON 排名不必与原文多模型排名一致，也不能据此说原文有误。

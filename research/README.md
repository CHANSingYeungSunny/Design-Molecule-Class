# 分子表示 Benchmark — 运行说明

本目录是本轮研究的全部成果。研究问题是：**固定同一个 LLM、同一批分子、同一道题和
同一套评分规则，只改变分子的输入表示，模型表现是否变化，Our Code 的优势和代价是什么？**

参考论文：*Rethinking Molecular Text Representations for LLMs: An Empirical Study*,
arXiv:2606.03057v1（2026-06-02）。本项目定位为
**「参考 MolRepBench 设计的扩展 benchmark」，不是论文复现**。

> **当前状态：框架已完成，smoke 规模真实调用已跑通（120/120 有效），main 规模尚未运行。**
> 主表里的数字来自 smoke（每格仅 8 题），**只能说明管线通畅，不能用于排名**。
> 正式比较需要先跑 main 规模，见 `07_benchmark/README.md`。

---

## 目录地图

| 目录 | 内容 | 状态 |
| --- | --- | --- |
| `00_baseline/` | 教授原版 `molecule_standalone.py`（只读基准记录）＋ `BASELINE.md` | 未修改 |
| `01_implementation/` | `molecule_core.py`：**交付核心，与基线逐字节相同** ＋ 撤销说明 | 本轮已还原 |
| `02_tests/` | `test_baseline_core.py`(55) ＋ `experiment_boundary.py`（边界探测） | 已运行通过 |
| `03_evidence/` | 源码证据摘录（E01–E23） | 历史材料 |
| `04_reference_notes/` | 参考项目抽取文本 | 历史材料 |
| `05_results/` | 边界探测输出、上一轮复现记录 | 历史材料 |
| `06_report/` | 旧报告（上一轮）＋ **本轮新汇报** | 见下 |
| `07_benchmark/` | **本轮新增：可运行的 benchmark 模块** | 已实现，待真实运行 |
| `08_deliverables/` | 项目实况报告、协议、manifest、进度矩阵、转换自检 | 本轮产出 |
| `99_archive_rollback_20261005/` | 核心撤销前的备份 | 备查 |

---

## 两个解释器，不要混用

| 用途 | 解释器 | 依赖 |
| --- | --- | --- |
| 教授原型核心 + 核心测试 | 系统 `python`（3.9.0），加 `-S -B` | **仅标准库** |
| benchmark 模块 + 报告 | `.venv/Scripts/python.exe`（3.13.5） | rdkit 2025.9.4、selfies、python-docx |

核心**只用标准库**是教授原型的硬约束，benchmark 侧用 RDKit 出参考标签是允许的
（RDKit 不进入核心类）。

---

## 一分钟版：不需要模型就能验证的完整流程

```powershell
$PY = "c:/Users/Asus/Desktop/Design Molecule Class/.venv/Scripts/python.exe"

# 1. 核心（教授原版）测试：55 项
cd "c:/Users/Asus/Desktop/Design Molecule Class/research/02_tests"
python -S -B -m unittest discover -s . -p "test_*.py"     # 预期 Ran 55 tests ... OK

# 2. benchmark 自检（不联网、不花钱）
cd ../07_benchmark
& $PY run_benchmark.py verify          # 三种表示往返覆盖率
& $PY run_benchmark.py selftest        # 评分器 17 项
& $PY run_benchmark.py offline-checks  # 客户端与密钥脱敏 9 项
& $PY run_benchmark.py dry-run         # 生成 120 题，不调用模型
```

本轮实测结果：55 项核心测试通过；三种表示往返覆盖率均 **1.000**；
评分器 **17/17**；离线检查 **9/9**；dry-run 生成 **120** 题。

---

## 需要模型才能跑的命令

```powershell
cd "c:/Users/Asus/Desktop/Design Molecule Class/research/07_benchmark"

# 先用 Our Code × 单任务 × 8 题打通链路
& $PY run_benchmark.py smoke --representations our_code --tasks t3_atom_counts --run-id probe-01

# 再放开三种表示
& $PY run_benchmark.py smoke --run-id smoke-01
& $PY run_benchmark.py main  --run-id main-01

# 从已保存响应重新评分 / 重新出表，不再花钱
& $PY run_benchmark.py score  --run main-01
& $PY run_benchmark.py tables --run main-01
```

现在直接跑会**拒绝启动**并列出缺少的配置项——这是设计行为。
需要配置的四项、地域与 base_url 对照表、以及本地配置方法，
见 `07_benchmark/README.md` 第 2 节。

---

## 本轮范围

**做**：Identity、官能团识别、重原子计数、**QED 估计（本项目新增任务）**、HBD 计数；
表示只做 **Our Code / SMILES / SELFIES** 三种。

**不做**（已登记但本轮不实现，注册表中标 `implemented: false`，不写空壳）：
图像、Graph2D、Graph3D、正式 MolJSON 对照；生成优化、docking、逆合成、模型训练。

**不使用 QM9**：不下载、不使用、不作为前置条件。数据用 **ChEBI-20**
（MolT5 再发布版，BSD-3-Clause）。

---

## 交付物

| 文件 | 用途 |
| --- | --- |
| `06_report/分子表示比较研究简明汇报.docx` | **本轮中文汇报**（面向教授） |
| `08_deliverables/build_report.py` | 汇报与 manifest 的生成入口 |
| `08_deliverables/project_audit_zh.md` / `.json` | 项目实况报告 |
| `08_deliverables/benchmark_protocol.md` | 数据、表示、任务、指标、公平性规则（动模型前冻结） |
| `08_deliverables/conversion_check.json` | 三种表示的往返一致性自检 |
| `08_deliverables/progress_matrix.csv` | 表示 × 任务的实现与运行状态 |
| `08_deliverables/experiment_manifest.json` | 代码、依赖、数据、模型、prompt、预算、种子、实际运行记录 |
| `07_benchmark/data/prepared/dataset_manifest.json` | 数据来源、校验值、筛选漏斗、标签支持度 |
| `07_benchmark/README.md` | benchmark 模块运行说明 |

---

## 已知边界（不要误读为支持）

- 核心不支持立体化学、CXSMILES 层、反应、查询语法、分子名称 → `UnsupportedSmilesError`
- 核心输出**确定性但非 canonical**
- 核心不做化学 sanitization（价键数、Kekulé、自由基不检查）
- 核心环闭合上限：同时打开的非树边 ≤ 99（原版行为，本轮未扩展）
- 构造代码往返：> `MAX_TEXT`(100000) 字符的图**不能**往返
- ChEBI-20 中 **68% 的分子带立体化学**，被排除而非剥掉，因此池只有 2,097 个分子；
  结论的适用范围是**无立体化学的小分子**
- 旧报告中的 16/16、4991 等数字**不是本轮结果**，对应材料不在工作区

## 历史材料

`03_evidence/`、`04_reference_notes/`、`05_results/`、`06_report/` 中上一轮的文档
保留为历史依据，**未篡改**。其中的数字不作为本轮证据引用。
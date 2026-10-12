# research/ 目录导航（当前状态）

> 本文件是**当前入口导航**。正式的运行步骤、结果表与限制写在根目录
> [README.md](../README.md)；审核结论在 [docs/AUDIT.md](../docs/AUDIT.md)；
> 论文对应在 [docs/PAPER_TO_PROJECT.md](../docs/PAPER_TO_PROJECT.md)；
> 核心代码索引在 [docs/CORE_CODE_INDEX.md](../docs/CORE_CODE_INDEX.md)。
>
> 本文件在 2026-10-13 之前描述的是**第 1 轮**（三种表示、`smoke` 阶段、`main` 未运行、
> 自检 17/17 与 9/9）。那份旧说明已原样归档到
> [docs/history/research_README_round1_2026.md](../docs/history/research_README_round1_2026.md)，
> 标注为历史阶段，**不代表当前状态**。

## 当前实验规模

| 项目 | 当前值 |
| --- | --- |
| 输入表示 | 8 种：SMILES、SELFIES、Our Code、img、Graph2D、Graph3D、CML、MolJSON |
| 任务 | 5 项：Identity、官能团识别、C/N/O/S 计数、QED 估计、HBD 计数 |
| 已完成运行 | `main-02`（原始日志 3000 行，全部 ok）、`extension-01`（原始日志 5038 行 = 5000 ok + 38 pending） |
| 当前 extension 冻结清单 | 4996 项 |
| 共同题分母 | Identity 200 题对，其余四项各 199 个分子；共同表逐题 7968 条 |
| 官能团 | 模型被问 **13 类**；主 Macro-F1 固定评分 **10 类**（aldehyde / nitrile / nitro 不进主分） |
| 模型 | 一个服务模型别名 `deepseek-v4.1-flash`；不训练模型、不调用工具 |

## 目录地图

| 目录 | 内容 | 说明 |
| --- | --- | --- |
| `00_baseline/` | 教授原版 `molecule_standalone.py` 与 `BASELINE.md` | 只读基准，字节未改（SHA-256 `55d314aa…`） |
| `01_implementation/` | `molecule_core.py`，与教授基线逐字节相同 | 生成 Our Code 文本；不参与评分 |
| `02_tests/` | `test_baseline_core.py`（55 项）、`experiment_boundary.py` | 核心测试，标准库即可运行 |
| `07_benchmark/` | benchmark 主体：`bmb/`、`config/`、`data/`、`runs/` | 见下"三条命令" |
| `08_deliverables/` | 共同题表生成器及其依赖、协议与 manifest | `build_common_table.py` 依赖 `build_analysis.py`、`build_extension_analysis.py` |
| `03_evidence/` `04_reference_notes/` `05_results/` `06_report/` `99_archive_*` | 早期材料与归档 | 历史材料，不作当前证据 |

## 三条命令各自验证什么

根目录 README 里的路径 A 有四条命令（`git clone`、`cd`、两条 `python`）。它们的作用分别是：

| 命令 | 作用 | 通过标志 |
| --- | --- | --- |
| `python research/02_tests/test_baseline_core.py` | 验证教授提供的 `Atom` / `Bond` / `Molecule` 基础类 | `Ran 55 tests ... OK` |
| `python scripts/replay_table.py` | **路径 A**：从保存的原始响应重解析答案，用原 `bmb/scoring.py` 重算 8×5 主表 | `Historical replay passed: 40 frozen lists, all 40 table values match.` |

**路径 B**（需要 RDKit）另有单独入口：

| 命令 | 作用 | 通过标志 |
| --- | --- | --- |
| `python scripts/rebuild_labels.py` | 用 RDKit 重新计算冻结题的参考答案（官能团、C/N/O/S、QED、HBD）并逐项比对 | `total 996 mismatched 0` |

只有路径 A 或只看教学例子 `scripts/demo_teaching.py`，**不能**说标签重建已通过。

## 历史内容在哪里

* 第 1 轮（三表示 / `smoke` 阶段）的说明：[docs/history/research_README_round1_2026.md](../docs/history/research_README_round1_2026.md)
* `07_benchmark/README.md` 正文同样以第 1 轮为主，文件开头已加当前状态提示，CLI 参数仍可参照。

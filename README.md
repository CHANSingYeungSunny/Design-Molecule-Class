# 分子表示比较实验：历史结果复算与复现

**研究问题**：同一个模型、同一批分子、同一批问题，只改变分子的**输入写法**，成绩会不会变？

我们用现有模型 `deepseek-v4.1-flash` 做了一次统一测试（benchmark）：把同一个分子写成 8 种输入形式，问 5 个化学问题，再用 RDKit（一个开源化学计算库）独立算出的参考答案打分。模型不调用工具、不执行代码、看不到参考答案。

**这是参考 MolRepBench 思路做的"扩展 benchmark"，不是那篇论文的完整复现。** 原论文是 9 种文本表示 × 8 项任务 × 16 个模型；本项目是 8 种输入形式 × 5 项任务 × 1 个模型，并且加上了图片输入、MOL 二维/三维文本、Our Code 和 QED 任务。逐条差异见 [docs/PAPER_TO_PROJECT.md](docs/PAPER_TO_PROJECT.md)。

**八种输入形式 × 五项任务**（↑ 越大越好，↓ 越小越好）：

| 表示 | 同分子判断 准确率 ↑ n=200 | 官能团 Macro-F1 ↑ n=199 | C/N/O/S 计数 全对率 ↑ n=199 | QED 估计 MAE ↓ n=199 | 氢键供体 准确率 ↑ n=199 |
| --- | ---: | ---: | ---: | ---: | ---: |
| SMILES | 0.985 | 0.839 | 0.231 | 0.111 | 0.814 |
| SELFIES | 0.540 | 0.546 | 0.015 | 0.180 | 0.432 |
| Our Code | 1.000 | 0.515 | 0.191 | 0.139 | 0.427 |
| img（二维结构图） | 0.990 | 0.868 | 0.342 | 0.119 | 0.975 |
| Graph2D（二维文本） | 0.850 | 0.438 | 0.055 | 0.145 | 0.382 |
| Graph3D（三维文本） | 0.980 | 0.571 | 0.070 | 0.170 | 0.387 |
| CML | 0.970 | 0.523 | 0.196 | 0.128 | 0.467 |
| MolJSON | 0.620 | 0.434 | 0.216 | 0.149 | 0.357 |

这八行是**八种写法**，不是八个 AI 模型。加粗/最好的数值只表示本列最高，**不表示已通过统计检验证明优于其他表示**。

---

## 三条复现路径

| 路径 | 做什么 | 需要什么 | 默认执行？ |
| --- | --- | --- | --- |
| **A 历史结果复算** | 从保存的模型原始回答重新解析、用原评分器重算 40 格 | Python 3.10+，不需要任何第三方包，不需要 API Key | ✅ 是 |
| **B 输入与标签重建** | 装完整环境，用 RDKit 重新计算冻结题的参考答案并逐项比对 | Python 3.13 + RDKit 等 | ✅ 是 |
| **C 重新调用模型** | 重新向服务发请求，产生新的回答 | 自备服务端点、API Key、预算 | ❌ 本任务不执行 |

---

## 路径 A：不用 API Key，重算整张表

```bash
git clone https://github.com/CHANSingYeungSunny/Design-Molecule-Class.git
cd Design-Molecule-Class

python research/02_tests/test_baseline_core.py
python scripts/replay_table.py
```

Windows 若没有 `python` 命令可用 `py -3`。

**成功标志**（按命令名称区分，不是按输出行数）

* `test_baseline_core.py`：验证教授提供的 `Atom` / `Bond` / `Molecule` 基础类，
  输出 `Ran 55 tests ... OK`。
* `replay_table.py`：路径 A 的历史复算，重新解析保存的模型原始回答并用原评分器重算 8×5 主表，
  输出 `Historical replay passed: 40 frozen lists, all 40 table values match.`
  并打印 `Denominators: {'t1_identity': 200, 't2_functional_groups': 199, 't3_atom_counts': 199, 't4_qed': 199, 't5_hydrogen_bond_donors': 199}`。

**输出目录 `reproduced/`**

| 文件 | 内容 |
| --- | --- |
| `summary.csv` | 重算出的 8×5 表格 |
| `per_item.jsonl` | 逐题证据：模型解析值、参考值、该题得分或误差 |
| `verification.json` | 校验结果：原始行数、排除行数、40 份冻结题单、图像字节数、覆盖率 |

这一步**不是**重新调用模型：它重新解析当时保存的模型原始 JSON，再调用原项目的 `bmb/scoring.py` 算分。它同时校验 40 份冻结题单的内容哈希、每条 prompt 的哈希、1200 张图片字节，以及八种表示在同一题上的对应关系。容差为 `1e-12`。

---

## 路径 B：完整环境，重建输入与参考答案

```bash
python -m venv .venv
# Windows PowerShell: .\.venv\Scripts\Activate.ps1
# macOS / Linux:      source .venv/bin/activate
python -m pip install -r requirements-benchmark.txt

python research/07_benchmark/run_benchmark.py selftest
python research/07_benchmark/run_benchmark.py offline-checks
python scripts/export_prompts.py
python scripts/rebuild_labels.py
python research/08_deliverables/build_common_table.py
python scripts/demo_teaching.py
```

**成功标志（本机 2026-10-13 实测，Windows CPython 3.13.5）**

| 检查 | 实测结果 |
| --- | --- |
| `research/02_tests/test_baseline_core.py` | 55 / 55 通过 |
| `run_benchmark.py selftest` | 99 checks, 0 failed |
| `run_benchmark.py offline-checks` | 14 offline checks, 0 failed |
| `scripts/check_code_index.py` | 53 行索引全部命中真实源码 |
| `scripts/export_prompts.py` | 导出 system prompt、8 种表示说明、5 个任务模板、5 个历史真实例子 |
| `scripts/rebuild_labels.py` | 996 条参考答案用 RDKit 重算，0 条不一致 |
| `scripts/replay_table.py` | 40 格全部一致 |
| `build_common_table.py` | 分母 200/199/199/199/199，完整性问题 0 |

**`scripts/rebuild_labels.py` 才是"标签重建"的验收**：它对 200 个 Identity 题对和 199×4 条理解题，用 RDKit 重新计算官能团（13 类 SMARTS 匹配）、C/N/O/S 计数、QED、HBD，并与冻结值逐项比对，同时检查八种表示在同一题上的分子与标签是否一致。只跑路径 A 或只看教学示例，**不能**宣称"完整标签重建通过"。

**教学示例**（明确不计入正式成绩）：`scripts/demo_teaching.py` 对乙醇 `CCO` 输出

```
Our Code : Molecule(atoms=(Atom(id='a0', element='C'), Atom(id='a1', element='C'),
                            Atom(id='a2', element='O')),
                    bonds=(Bond(begin='a0', end='a1'), Bond(begin='a1', end='a2')))
SMILES   : CCO
SELFIES  : [C][C][O]
RDKit 参考: C=2, N=0, O=1, S=0；HBD=1；QED=0.406808
```

这个分子太小（不满足"分子量 ≥ 300 Da 且至少 2 个环"的门槛），**不进入正式成绩**。任何地方都不使用伪造的模型答案冒充历史成绩。

**版本记录**（实测，见 `requirements-benchmark.txt`）

* Python 3.13.5
* rdkit：安装包版本 `2025.9.4`，`rdkit.__version__` 报告 `2025.09.4`
* selfies：安装包版本 `2.2.0`，但 `selfies.__version__` 报告 `2.1.1`（该分发包升级时没有同步模块内的版本号，两个数都记录，避免误判）
* numpy 2.5.3、Pillow 12.3.0、python-docx 1.2.0

**不要对历史快照运行 `prepare`、`main-prep`、`extension-prep`**：这些命令会重新生成冻结文件。探索新设计请用独立分支和新的 run ID。

---

## 路径 C：重新调用模型（本任务默认不执行）

* 需要自备适用的服务端点、API Key 与预算。历史 run 记录为：`provider=token_plan`、端点 `https://token-plan.cn-beijing.maas.aliyuncs.com/compatible-mode/v1`、模型 `deepseek-v4.1-flash`、`temperature=0`、`top_p=1`、`max_tokens=4096`、zero-shot、`extra_body.enable_thinking=false`。
* 密钥只放环境变量或本地被忽略的 `research/07_benchmark/config/local.env`；`local.env.example` 只有变量名。可用性与费用限制由真实服务配置决定，本仓库不代付、不代申请。
* 真实可用的入口是 `research/07_benchmark/run_benchmark.py`（用 `python research/07_benchmark/run_benchmark.py --help` 查看实际支持的子命令与参数），请另起新的 run ID。
* 即使 `temperature=0`，服务端权重版本或推理实现变化都可能让新回答不同，**不能承诺逐题一致**。只有路径 A 的 replay 能精确核对旧表。新模型/新端点/新 prompt 的结果必须独立报告。

---

## 文件导航

| 想看什么 | 去哪里 |
| --- | --- |
| 研究背景、三条路径、结果表 | 本 `README.md` |
| 原论文逐条对应与差异 | [docs/PAPER_TO_PROJECT.md](docs/PAPER_TO_PROJECT.md) |
| 核心函数真实路径、行号、中文说明 | [docs/CORE_CODE_INDEX.md](docs/CORE_CODE_INDEX.md) |
| 审核结论、发现的限制、验证状态 | [docs/AUDIT.md](docs/AUDIT.md) |
| 向教授口头说明与追问清单 | [docs/EXPLAIN_TO_PROF.md](docs/EXPLAIN_TO_PROF.md) |
| 被测模型的 system prompt | [prompts/system.txt](prompts/system.txt) |
| 八种表示说明 + 五项任务模板 | [prompts/benchmark_templates.json](prompts/benchmark_templates.json) |
| 五项任务的真实历史请求与回答 | [prompts/historical_our_code_examples.json](prompts/historical_our_code_examples.json) |
| 基础类 Atom / Bond / Molecule | [research/01_implementation/molecule_core.py](research/01_implementation/molecule_core.py) |
| 实际模型 prompt 组装 | [research/07_benchmark/bmb/tasks.py](research/07_benchmark/bmb/tasks.py) |
| 独立 RDKit 参考答案 | [research/07_benchmark/bmb/labels.py](research/07_benchmark/bmb/labels.py) |
| 评分规则 | [research/07_benchmark/bmb/scoring.py](research/07_benchmark/bmb/scoring.py) |
| 无 Key 历史复算入口 | [scripts/replay_table.py](scripts/replay_table.py) |
| 参考答案重建入口 | [scripts/rebuild_labels.py](scripts/rebuild_labels.py) |
| 原始模型回答（历史证据） | `research/07_benchmark/runs/main-02/raw_responses.jsonl`、`research/07_benchmark/runs/extension-01/raw_responses.jsonl` |

**两种 prompt 不要混淆**：`prompts/` 下的是**被测模型**实际使用的模板，直接导出自 `tasks.py`。`docs/codebuddy_execution_prompt.txt` 是**开发工具 CodeBuddy** 的工作指令，不是分子题目。

**历史格式缺陷已保留**：T2/T3 的真实 prompt 示例中含双大括号（例如 `{{"C": 0}}`）。这是当时就存在的格式缺陷，历史 prompt 与响应一律不改。若未来修正，必须另立 prompt 版本、在 dev 上验证、使用新的 run ID，不得改写旧记录。

---

## 已知限制

1. **Identity 的 1.000 只作描述性成绩。** 原三表示的正例用视图 0/1、负例用 0/0，新增五表示两类都用 0/1；历史输入排序可能带标签线索，题对还可能共享分子端点。不能据此声称 Our Code 的同分子判断能力已确认优于其他表示。
2. **13 类被问、10 类被评分。** 官能团实际问 13 类，主 Macro-F1 只对冻结的 10 类平均；aldehyde、nitrile、nitro 被问到但不进主分。
3. **原始行数 ≠ 有效题数。** `main-02` 原始日志 3000 行（全部 `ok`）；`extension-01` 原始日志 5038 行，其中 38 条是 `pending` 占位、4 条属于已被当前冻结清单取代的历史题。当前冻结清单为 4996 项。复算按当前清单过滤、按 `example_id` 取最后状态。
4. **少 1 个分子有明确原因。** 理解任务分母是 199 而不是 200：一个含铜分子不满足当前 Graph3D 构象流程。这是转换覆盖问题，必须另报，不能当成模型答错，也不能悄悄换一题。
5. **QED / HBD 是结构可计算指标**，不能推断疗效或安全性。
6. **Our Code + H 只有结构一致性验证**，尚未进行模型消融实验，没有新成绩。
7. **重新推理不一定逐题一致**，见路径 C。
8. **只有一个服务模型别名**，服务端未提供权重版本；缺失 `reasoning` 字段时思考状态仍不可确认。

---

## 出处与许可

* 教授提供的基础类（`Atom` / `Bond` / `Molecule`）保持原字节，三份副本 SHA-256 均为 `55d314aa312707b498536d339ff59dc89bbbe4956cab54a74a4746cd385af57b`，来源为 `research/00_baseline/`。本报告与代码不把它说成新的原创实现。
* `reference/MolJSON/` 保留第三方 MolJSON 的源码、schema 与其 BSD-3-Clause LICENSE（来源 oxpig/MolJSON，上游 commit 记录见 `docs/AUDIT.md`）。不为第三方内容改换许可证。
* RDKit（BSD-3-Clause）通过 pip 安装，不在本仓库内附带整份源码。
* 分子数据来自 ChEBI-20 的 MolT5 再发布文件（`research/07_benchmark/data/raw/`），划分由本项目重做。
* 本仓库没有为整体项目选定新的开源许可证；第三方内容各自保留原许可。

---

## 发布状态（实测）

| 项目 | 实际值 |
| --- | --- |
| 仓库 | <https://github.com/CHANSingYeungSunny/Design-Molecule-Class> |
| **已验证源码提交** | `aac92e78a7d1b1f80653d59ec45252a3a612e45a` |
| 验证方式 | 从全新目录 clone 该 commit：`test_baseline_core.py` → `Ran 55 tests ... OK`；`replay_table.py` → `Historical replay passed: 40 frozen lists, all 40 table values match.`；`rebuild_labels.py` → 996 条参考答案、0 条不一致；`reference/MolJSON` 源码在 clone 后可正常 import |
| 远端自动化检查 | GitHub Actions 工作流 `Historical offline replay`（ubuntu-24.04，Python 3.13，无任何 secret）。已核对 `head_sha` 的两条成功运行分别见下表 |
| 文档与报告提交 | 见仓库 Commits 页面最新一条；只新增/更新文档、脚本与报告，不改动核心代码、冻结题单与原始响应 |
| 未执行 | 重新付费调用模型（路径 C）；需要自备端点、key 与预算 |

**已核对 `head_sha` 的 Actions 运行**（用 `gh api .../actions/runs/<id>` 读取，不是靠"有两个成功 run"推断）

| run ID | 对应 head_sha | 状态 | 结论 | 链接 |
| --- | --- | --- | --- | --- |
| 38315294795 | `aac92e78a7d1b1f80653d59ec45252a3a612e45a` | completed | success | <https://github.com/CHANSingYeungSunny/Design-Molecule-Class/actions/runs/38315294795> |
| 38316032354 | `10ef071c6f46c57f1118fe085a72e7fc55c23687` | completed | success | <https://github.com/CHANSingYeungSunny/Design-Molecule-Class/actions/runs/38316032354> |

最终交付提交的运行会在推送后另行核对，不能用上表代替本次验收。

**关于"谁做的验收"**：上面的 clone 复算是在**与原工作目录分离的全新目录**由本次工作完成的，
**没有另一位独立使用者的验收记录**。仓库提供了 README 与命令，可供其他人自行复算；
但"可供他人复算"与"已经由其他人独立验证"是两件事，本文件与报告都按这个界线表述。

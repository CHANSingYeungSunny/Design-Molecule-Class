# 审核结论、已验证项与限制

**验证日期：2026-10-13（香港时间）。本文件只写本机与远端实际执行过的检查，未执行的一律标注"未执行"。**

## 1 仓库与发布状态

| 项目 | 实际值 |
| --- | --- |
| 目标仓库 | https://github.com/CHANSingYeungSunny/Design-Molecule-Class |
| 已验证源码提交 | `aac92e78a7d1b1f80653d59ec45252a3a612e45a`（fresh clone 与 Actions 都针对这个 commit 执行） |
| 文档与本报告提交 | 见仓库 Commits 页面最新一条；该提交只新增/更新文档、脚本与报告，不改动核心代码、冻结题单与原始响应 |
| 远端初始状态 | 推送前用 `gh api` 确认为**空仓库**（`/commits` 409 "Git Repository is empty"），因此首次推送不会覆盖任何已有工作；未使用 force push |
| 分支 | 只使用 `main`；没有 force push、没有 `reset --hard`、没有改写历史 |

**已用 `gh api .../actions/runs/<id>` 核对 `head_sha` 的运行**（不靠"有两个成功 run"推断对应关系）

| run ID | 对应 head_sha | status | conclusion | 链接 |
| --- | --- | --- | --- | --- |
| 38315294795 | `aac92e78a7d1b1f80653d59ec45252a3a612e45a` | completed | success | <https://github.com/CHANSingYeungSunny/Design-Molecule-Class/actions/runs/38315294795> |
| 38316032354 | `10ef071c6f46c57f1118fe085a72e7fc55c23687` | completed | success | <https://github.com/CHANSingYeungSunny/Design-Molecule-Class/actions/runs/38316032354> |

最终交付提交的运行会在推送后另行核对；不能用上表代替最终验收。

**"谁做的验收"的界线**：下面的 clone 复算是在**与原工作目录分离的全新目录**完成的，
执行者仍是本次工作，**没有另一位独立使用者的验收记录**。仓库提供 README 与命令可供其他人自行复算，
但"可供他人复算"不等于"已经由其他人独立验证"。

## 2 本机实测（Windows，`\.venv\Scripts\python.exe`，CPython 3.13.5）

| 检查 | 命令 | 实测结果 |
| --- | --- | --- |
| 教授基础类 | `research/02_tests/test_baseline_core.py` | `Ran 55 tests ... OK` |
| benchmark 自检 | `research/07_benchmark/run_benchmark.py selftest` | `99 checks, 0 failed` |
| 客户端离线检查 | `research/07_benchmark/run_benchmark.py offline-checks` | `14 offline checks, 0 failed` |
| 核心代码索引 | `scripts/check_code_index.py` | `index rows checked: 53` / `core code index OK` |
| prompt 导出 | `scripts/export_prompts.py` | system prompt、8 种表示说明、5 个任务模板、5 个真实历史例子 |
| **路径 A：历史复算** | `scripts/replay_table.py` | `Historical replay passed: 40 frozen lists, all 40 table values match.`；分母 200/199/199/199/199 |
| 项目自带的共同题表入口 | `research/08_deliverables/build_common_table.py` | 分母 200/199/199/199/199，完整性问题 0，40 格数值与 replay 一致 |
| **路径 B：参考答案重建** | `scripts/rebuild_labels.py` | RDKit 2025.09.4 重算 996 条（200 题对 + 199×4），0 条不一致；同时校验八臂分子与标签一致 |
| 教学示例 | `scripts/demo_teaching.py` | CCO → Our Code / `CCO` / `[C][C][O]`；C=2,N=0,O=1,S=0；HBD=1；QED=0.406808（明确不计入正式成绩） |

原始日志的行数分布（逐行统计，不是引用旧数字）：

* `main-02/raw_responses.jsonl`：3000 行，全部 `ok`。
* `extension-01/raw_responses.jsonl`：5038 行 = 5000 `ok` + 38 `pending`；另有 4 条属于已被当前冻结清单取代的历史 Identity 题（复算按当前清单过滤时自动排除）。

`build_common_table.py` 会重新写出 `research/08_deliverables/analysis_common_items.json`。本机在 Windows 上重跑后文件内容与冻结版本**语义完全一致**，但换行被写成 CRLF；已把冻结版本的字节恢复原样，避免历史交付物被换行改写。

## 3 路径 A / B / C 的界线（不要混淆）

* **A（默认、已完成、不付费）**：从保存的原始响应重解析答案，用原 `bmb/scoring.py` 重算 40 格。它验证"表格可以从保存数据重建"。
* **B（已完成、不付费）**：用 RDKit 对冻结题目**重新计算参考答案**并逐项比对。它验证"参考答案可复现"。只有 A 或只有教学例子，**不能**宣称标签重建通过。
* **C（未执行）**：重新向模型服务发请求。需要自备端点、key、预算；即使 `temperature=0`，服务端权重或推理实现变化都可能让新回答不同，不能承诺逐题一致。本任务未获授权，也没有执行。

## 4 与 ChatGPT 候选交付包的关系

用 SHA-256 逐文件比对候选包与本机：19 个核心文件里 17 个逐字节一致；仅两处不同，均已审阅并采纳：

1. `bmb/runner.py`：新增 `import pathlib`，把冻结图片相对路径改为 `PureWindowsPath(...).parts` 解析，使 Windows 冻结的路径在 Linux 也能读取（改动仅此两行，未触碰历史数据）。
2. `run_benchmark.py`：离线脱敏测试的假密钥字符串改成明确的 `sk-offline-selftest-dummy-000000`，不用于任何服务请求。

其余本机更晚或更完整的内容一律保留本机版本，例如 `research/README.md`（本地 6241 字节，比候选包的 404 字节完整）。候选包的 prompt 导出文件与本机重新导出的文件**逐字节一致**（修正了 `export_prompts.py` 在 Windows 上把换行写成 CRLF 的问题后）。

## 5 prompt 洁净性（代码路径证据）

* `bmb/tasks.py` 的 `system_message()`（137–143 行）与五个任务的 `messages(rep, item, run_config)`（275 / 333 / 398 / 461 / 524 行）函数体只读 `item['payload']` 与常量说明文本，**从不读取 `item['label']`、参考答案、CID、分子名称**；`example_id` 不进入消息正文；`run_config` 未参与拼装。
* 参考答案只出现在 `runner.py` 第 477 / 484 / 493 / 513 / 528 行附近写入记录并送进评分，与实际发送的 `messages` 在代码路径上分离。
* img 通过 `tasks.py` 72–89 行的 `image_url` 内容块（`data:image/png;base64,...`）发送，不是文字描述；PNG 经 Pillow 重编码后不含文本元数据块（自检断言见 `run_benchmark.py` 1130–1157 行）。
* 请求体 `bmb/client.py` `_payload`（254–266 行）只含 `model / messages / temperature / top_p / max_tokens / stream / response_format`，**没有 `tools`**（自检断言见 `run_benchmark.py` 686 行）。
* 全仓库没有 `eval` / `exec` / `subprocess` 执行 Our Code 文本的通路。`core_loader.py` 用 `importlib` 加载基础类只是让实验程序自己生成表示与做往返校验；被测模型一侧没有任何工具、代码执行或检索通道。
* 需要披露的一点：`FunctionalGroupTask.messages` 列出的基团名与答案的 JSON 键同名——这是答题**格式（键名）**，不是答案**值（true/false）**。

## 6 已知限制（保留并披露）

1. **Identity 只作描述性比较。** 原三表示正例视图 0/1、负例 0/0；新增五表示两类均 0/1。历史输入排序可能带标签线索，题对可能共享分子端点（`tasks.py` 247–254 行）。不能据此说 Our Code 的同分子判断能力已确认优于其他表示。
2. **13 类被问、10 类被评分。** aldehyde / nitrile / nitro 被问到但不进主 Macro-F1。
3. **T2/T3 历史示例含双大括号**（例如 `{{"C": 0}}`）。历史 prompt 与响应一律不改；未来修正必须另立 prompt 版本、在 dev 验证、使用新的 run ID。
4. **历史通道记录**：run manifest 记录 `provider=token_plan`、端点 `https://token-plan.cn-beijing.maas.aliyuncs.com/compatible-mode/v1`、模型 `deepseek-v4.1-flash`、`temperature=0`、`top_p=1`、`max_tokens=4096`、zero-shot、`extra_body.enable_thinking=false`（`main-02.log` 前 6 行可直接核对）。早期笔记里的"普通 DashScope 端点"说法不适用于这次历史运行。注意 `config/local.env.example` 描述的是普通 API 的申请方式，与服务端条款相关；缺失 `reasoning` 字段时思考状态仍不可确认，服务端也没有固定权重版本的证据。
5. **HBD 口径**：prompt 用"带氢的 N/O 位点"作简化解释，oracle 用 RDKit `CalcNumHBD`；特殊原子的定义边界应在新版本检查，不据此改历史标签。
6. **转换覆盖**：理解任务分母是 199 而非 200，因为一个含铜分子不满足当前 Graph3D 构象流程。这是转换问题，必须另报，不能当模型答错，也不能换题。
7. **QED / HBD 是结构可计算指标**，不能推断疗效或安全性。
8. **Our Code + H 只有结构一致性验证**，没有模型消融成绩。

## 7 发布卫生

* 未提交：`.venv/`、`research/07_benchmark/config/local.env`（真实密钥）、`reference/rdkit/`（约 102 MiB 第三方整仓）、`extension_items_v1/`（约 64 MiB 重复）、非发布 run 目录、`scored.jsonl`、`companion.csv`、`*.bak`、临时比对文件。
* 最大单文件为 `runs/extension-01/raw_responses.jsonl` = 89.7 MiB（94,061,418 字节），低于 GitHub 单文件上限，因此**未使用 Git LFS**；本文件不写入任何未核实的 LFS 额度说明。
* 已提交前扫描：暂存 1766 个文件中无任何形如密钥的字符串；`local.env.example` 只含变量名。
* 换行策略：`.gitattributes` 对 `*.jsonl`、`*.png`、`*.docx`、核心 `.py` 与 `prompts/` 设 `-text`，防止换行改写导致字节哈希失配；其余文本统一 LF。
* 第三方：`reference/MolJSON` 保留其 BSD-3-Clause LICENSE，上游为 oxpig/MolJSON，commit `047b135e628b091ec57a65d6da06989554f3a229`（读取其本地 `.git` 记录，未改动该目录）。`moljson` 包被 `bmb/reps_ext.py` 的 `moljson_text` 真实 import，因此源码必须随仓库分发。
* 本仓库没有为整体项目选定新的开源许可证；第三方内容各自保留原许可。

## 8 报告渲染与页数

`分子表示比较研究简明汇报_过程与复现版_GitHub验证版.docx` 在**本机 Microsoft Word（Office 16）导出 PDF 时为 9 页**；
同一份 DOCX 在另一个渲染环境下曾得到 10 页。DOCX 是流式排版格式，分页取决于渲染器的字体度量与页面设置，
**页数不是验收标准**：验收看的是章节不重复、40 格成绩未改、分母 200/199 与 13 问 10 评未变、链接可查。
仓库里同时提供 DOCX 与用本机 Word 导出的 PDF；PDF 只是为了固定分页便于阅读与打印。

文中展示的源码片段为排版做了换行处理，**打印换行不等于可执行的完整源码**；完整源码以仓库中对应文件的
固定链接为准。

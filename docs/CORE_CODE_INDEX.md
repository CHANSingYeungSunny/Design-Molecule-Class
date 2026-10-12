# 核心代码位置索引

行号对应当前仓库源码，可由 `python scripts/check_code_index.py` 自动复核（它会把下表每一行的路径与行号读回源码验证，44 行全部通过）。
正式发布后应绑定实际 commit：行号只对某一版本有效，固定链接才是永久证据。

## 1 行号索引

| 文件 | 定义 | 起止行 |
| --- | --- | --- |
| `research/01_implementation/molecule_core.py` | `Atom` | 43–84 |
| `research/01_implementation/molecule_core.py` | `Bond` | 88–97 |
| `research/01_implementation/molecule_core.py` | `Molecule` | 351–489 |
| `research/01_implementation/molecule_core.py` | `Molecule.to_code` | 438–449 |
| `research/01_implementation/molecule_core.py` | `Molecule.from_code` | 452–489 |
| `research/07_benchmark/bmb/dataset.py` | `read_source_rows` | 49–60 |
| `research/07_benchmark/bmb/dataset.py` | `build_pool` | 63–153 |
| `research/07_benchmark/bmb/dataset.py` | `split_pool` | 156–166 |
| `research/07_benchmark/bmb/dataset.py` | `apply_complexity_filter` | 275–290 |
| `research/07_benchmark/bmb/reps.py` | `canonical_smiles` | 47–56 |
| `research/07_benchmark/bmb/reps.py` | `smiles_text` | 69–75 |
| `research/07_benchmark/bmb/reps.py` | `selfies_text` | 95–106 |
| `research/07_benchmark/bmb/reps.py` | `our_code_text` | 122–127 |
| `research/07_benchmark/bmb/reps.py` | `our_code_identity_variant` | 130–149 |
| `research/07_benchmark/bmb/reps.py` | `render` | 191–206 |
| `research/07_benchmark/bmb/reps_ext.py` | `img_png` | 205–207 |
| `research/07_benchmark/bmb/reps_ext.py` | `graph2d_text` | 351–362 |
| `research/07_benchmark/bmb/reps_ext.py` | `graph3d_text` | 506–553 |
| `research/07_benchmark/bmb/reps_ext.py` | `cml_text` | 579–580 |
| `research/07_benchmark/bmb/reps_ext.py` | `moljson_text` | 660–662 |
| `research/07_benchmark/bmb/tasks.py` | `system_message` | 137–143 |
| `research/07_benchmark/bmb/tasks.py` | `user_message` | 77–89 |
| `research/07_benchmark/bmb/tasks.py` | `payload_text` | 92–96 |
| `research/07_benchmark/bmb/tasks.py` | `IdentityTask.build_pairs` | 170–235 |
| `research/07_benchmark/bmb/tasks.py` | `IdentityTask.messages` | 275–283 |
| `research/07_benchmark/bmb/tasks.py` | `FunctionalGroupTask.messages` | 333–346 |
| `research/07_benchmark/bmb/tasks.py` | `AtomCountTask.messages` | 398–410 |
| `research/07_benchmark/bmb/tasks.py` | `QedTask.messages` | 461–472 |
| `research/07_benchmark/bmb/tasks.py` | `DescriptorTask.messages` | 524–534 |
| `research/07_benchmark/bmb/tasks.py` | `prompt_hash` | 569–570 |
| `research/07_benchmark/bmb/labels.py` | `functional_group_labels` | 141–160 |
| `research/07_benchmark/bmb/labels.py` | `atom_count` | 163–167 |
| `research/07_benchmark/bmb/labels.py` | `atom_counts` | 170–171 |
| `research/07_benchmark/bmb/labels.py` | `qed_reference` | 174–178 |
| `research/07_benchmark/bmb/labels.py` | `descriptor_reference` | 181–193 |
| `research/07_benchmark/bmb/labels.py` | `is_same_molecule` | 196–208 |
| `research/07_benchmark/bmb/labels.py` | `frozen_macro_f1_labels` | 102–122 |
| `research/07_benchmark/bmb/client.py` | `ModelClient._payload` | 254–266 |
| `research/07_benchmark/bmb/client.py` | `ModelClient.chat` | 268–316 |
| `research/07_benchmark/bmb/runner.py` | `load_frozen_items` | 238–295 |
| `research/07_benchmark/bmb/runner.py` | `load_extension_items` | 298–356 |
| `research/07_benchmark/bmb/runner.py` | `run_cell` | 382–537 |
| `research/07_benchmark/bmb/runner.py` | `score_run` | 615–662 |
| `research/07_benchmark/bmb/scoring.py` | `score_item` | 86–101 |
| `research/07_benchmark/bmb/scoring.py` | `_functional_group_metrics` | 137–220 |
| `research/07_benchmark/bmb/scoring.py` | `_atom_count_metrics` | 223–243 |
| `research/07_benchmark/bmb/scoring.py` | `_qed_metrics` | 246–265 |
| `research/07_benchmark/bmb/scoring.py` | `_descriptor_metrics` | 268–279 |
| `research/07_benchmark/bmb/scoring.py` | `pair_key` | 649–666 |
| `research/08_deliverables/build_common_table.py` | `main` | 54–196 |
| `scripts/replay_table.py` | `main` | 63–200 |
| `scripts/rebuild_labels.py` | `main` | 103–173 |
| `scripts/export_prompts.py` | `main` | 20–54 |

## 2 每个环节在做什么

### 基础类：`research/01_implementation/molecule_core.py`

* 干什么：用 `Atom`（一个原子）/ `Bond`（两个原子之间的键）/ `Molecule`（把它们收成一个分子）三个记录类描述分子，`Molecule.to_code()` 把这个分子打印成一段构造代码文本，也就是被测模型看到的 "Our Code"。
* 输入：`to_code()` 不需要外部输入，它把已经建好的 `Molecule` 对象输出成字符串；`from_code()` 反向读字符串。
* 输出：形如 `Molecule(atoms=(Atom(id='a0', element='C'), ...), bonds=(Bond(begin='a0', end='a1'), ...))` 的文本。
* 注意：这三份文件（根目录 `molecule_standalone.py`、`research/00_baseline/molecule_standalone.py`、`research/01_implementation/molecule_core.py`）字节完全一致，SHA-256 都是 `55d314aa312707b498536d339ff59dc89bbbe4956cab54a74a4746cd385af57b`，来自教授提供的原型，本次未改动。它用来**生成**表示，不参与评分——评分一律走 RDKit。

### 数据筛选：`research/07_benchmark/bmb/dataset.py`

* 干什么：读取 ChEBI-20 的再发布 validation/test 文件，合并后重划 dev/test，并施加复杂度门槛（分子量 ≥ 300 Da、环数 ≥ 2）与 Our Code 可解析/往返一致的限制。
* 输入：原始 txt 数据文件；输出：候选池条目（SMILES、来源、划分）。
* 当前复杂度候选池 431 个分子。dev 只用于调方案，不进入正式成绩。

### 表示转换：`bmb/reps.py` 与 `bmb/reps_ext.py`

* 干什么：把同一个参考分子分别写成 8 种输入形式。
* 输入：一个规范 SMILES；输出：该表示下的文本，或（img）PNG 字节。
* 分工：`reps.py` 负责 SMILES / SELFIES / Our Code；`reps_ext.py` 负责 img（RDKit 绘制后经 Pillow 重编码，剥掉 PNG 文本元数据）、Graph2D / Graph3D（MOL 连接表文本，3D 坐标来自计算构象）、CML（XML）、MolJSON（JSON，来自第三方 `moljson` 包）。
* 教学例子（不计入正式成绩）：乙醇 `CCO` → Our Code 为 3 个 `Atom` + 2 个 `Bond`，SMILES `CCO`，SELFIES `[C][C][O]`。

### 模型 prompt：`bmb/tasks.py`

* 干什么：组装真正发给被测模型的 system / user 消息。
* 输入：`rep`（表示名）、`item['payload']`（该表示下的分子内容）、`run_config`（实际未参与拼装）。
* 输出：`[system_message(), user_message(...)]`，即两到三条消息。
* 已核对的洁净性：五个 `messages()` 的函数体只读 `item['payload']` 与常量说明文本，**从不读 `item['label']`、参考答案、CID 或分子名称**；`example_id` 也不进入消息正文。实际发出的提示词原文会随响应一起落盘为 `prompt_text`，可逐条审计。
* img 走 `tasks.py:72–89` 的 `image_url` 内容块（`data:image/png;base64,...`），不是文字描述。
* Our Code 只是给模型**阅读**的文本。全仓库没有 `eval` / `exec` / `subprocess` 去执行它；`core_loader.py` 用 `importlib` 加载基础类只是为了让实验程序自己生成表示、做往返校验，与模型无关。

### 独立参考答案（oracle）：`bmb/labels.py`

* 干什么：用 RDKit 独立算出参考答案，也就是"标准答案"。
* 输入：参考分子的规范 SMILES；输出：官能团 13 类布尔值、C/N/O/S 计数、QED 浮点值、HBD 整数、Identity 的同分子判定。
* 关键：它**不使用**我们的基础类，所以 Our Code 不能给自己的表示打分。
* 本仓库新增 `scripts/rebuild_labels.py` 会用 RDKit 把这 996 条参考答案全部重算一遍并与冻结值逐项比对。

### 请求与日志：`bmb/client.py` 与 `bmb/runner.py`

* `client.py` 的 `_payload`(254–266) 构造请求体：只含 `model / messages / temperature / top_p / max_tokens / stream / response_format`，**没有 `tools`**，也没有任何函数调用通道。
* `runner.py` 的 `run_cell`(382–537) 读冻结题单、逐条调用、把原始回答与 `prompt_text`、`label`、`usage`、`prompt_hash` 一起追加写入 JSONL，可断点续跑。
* `label` 只在这几处被写入记录并送进评分，从不进入 `messages`。

### 评分：`bmb/scoring.py`

* 干什么：把模型的解析结果和参考答案比成分数。
* `score_item`(86–101) 是逐题判定；`_functional_group_metrics`(137–220) 先按每一类统计正确识别/误报/漏报，算 F1 再等权平均；`_atom_count_metrics`(223–243) 要求 C/N/O/S 四个数全对；`_qed_metrics`(246–265) 取绝对误差的均值（MAE）；`_descriptor_metrics`(268–279) 用于 HBD 的精确匹配。
* 输入：一批响应记录 + 冻结标签集；输出：每个表示×任务的指标值。

### 共同题表：`research/08_deliverables/build_common_table.py`

* 干什么：对每项任务取 8 种表示的共同题交集，然后重新算指标，得到 8×5=40 格。
* 输入：`runs/main-02`、`runs/extension-01` 的原始响应与冻结清单（依赖 `build_analysis.py` 与 `build_extension_analysis.py`，不要只提交这一个脚本）。
* 输出：`research/08_deliverables/analysis_common_items.json`。
* 另一个入口 `scripts/replay_table.py` 用纯标准库做同样的事，无 API Key 也能跑。

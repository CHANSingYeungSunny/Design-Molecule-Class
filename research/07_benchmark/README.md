# 07_benchmark — 分子表示对比评测模块

**参考 MolRepBench 设计的扩展 benchmark**，不是原文复现。参考论文：
*Rethinking Molecular Text Representations for LLMs: An Empirical Study*,
arXiv:2606.03057v1（2026-06-02），作者 Arun Raja, Garrett M. Morris, Kian Ming A. Chai。

> **当前状态提示（2026-10-13 加）**：本文件正文记录的是**第 1 轮**（三种表示、
> `smoke` 阶段、`main` 未运行、自检 17/17 与 9/9）。后续阶段已经把表示扩到
> **八种**、任务扩到 **五项**，并且已经跑完 `main-02` 与 `extension-01`。
> 下面的"第 1 轮"描述与 `smoke/probe` 相关状态**是历史记录，不是当前状态**；
> 当前结果、复现命令与限制请看根目录 [README.md](../../README.md) 与
> [docs/AUDIT.md](../../docs/AUDIT.md)。CLI 参数与模块结构仍以本文件为准。
>
> 第 1 轮的目录级说明已归档到
> [docs/history/research_README_round1_2026.md](../../docs/history/research_README_round1_2026.md)。

本轮（第 1 轮）实测三种表示：**Our Code、SMILES、SELFIES**。图像、Graph2D、
Graph3D、MolJSON 已在注册表中登记但**本轮不实现**，注册表里标 `implemented: false`
并写明原因，不写空壳函数。

---

## 0. 环境与解释器（重要）

本项目有两个解释器，用途不同，不要混用：

| 用途 | 解释器 | 依赖 |
| --- | --- | --- |
| 教授原型核心（`01_implementation/molecule_core.py`） | 系统 `python`（3.9.0） | **仅标准库** |
| 本 benchmark 模块 | `.venv/Scripts/python.exe`（3.13.5） | rdkit 2025.9.4、selfies、python-docx |

核心是教授原版，逐字节等于 `00_baseline/molecule_standalone.py`（SHA-256
`55d314aa…af57b`）。benchmark **不修改**核心，也不通过核心算分：参考标签全部由
RDKit 产生，核心只负责把分子渲染成模型要读的 `Molecule(...)` 文本。

以下命令一律在 `research/07_benchmark/` 目录下执行，`PY` 指 venv 解释器：

```powershell
cd "c:/Users/Asus/Desktop/Design Molecule Class/research/07_benchmark"
$PY = "c:/Users/Asus/Desktop/Design Molecule Class/.venv/Scripts/python.exe"
```

---

## 1. 不需要模型就能跑的命令（本轮已实测）

```powershell
& $PY run_benchmark.py prepare        # 下载数据已就位时：建样本池、划分、写 manifest
& $PY run_benchmark.py verify         # 三种表示对全池的往返一致性 + Identity 题对自检
& $PY run_benchmark.py selftest       # 评分器自检（已知答案 + 故意错误输出）
& $PY run_benchmark.py offline-checks # 客户端请求形状与密钥脱敏检查，不联网
& $PY run_benchmark.py dry-run        # 生成全部 120 道题与 prompt，不调用模型
```

本轮实际结果见 `../08_deliverables/`。要点：

- `prepare`：6,601 行 → 池 **2,097** 个分子（dev 524 / test 1,573）。
- `verify`：三种表示往返覆盖率均 **1.000**；Identity 正例 32/32 同分子、负例 32/32 异分子。
- `selftest`：**17 项全过**。
- `offline-checks`：**9 项全过**。

---

## 2. 需要模型配置才能跑的命令

```powershell
& $PY run_benchmark.py smoke          # 真实调用，每任务 8 题
& $PY run_benchmark.py main           # 真实调用，每任务 200 题
& $PY run_benchmark.py score --run <run_id>     # 从已保存响应重新评分，不再花钱
& $PY run_benchmark.py tables --run <run_id>    # 从结果文件重新生成表格
```

**当前状态**：`smoke` 已跑通（run `smoke-02`，120/120 有效，0 截断）；
**`main` 未运行**。未配置端点时 `smoke`/`main` 会**拒绝启动**并列出缺少的配置项——这是设计行为。

### 2.0 实验使用的端点（2026-10-06 迁移后）

| 项目 | 值 |
| --- | --- |
| base_url | **`https://dashscope.aliyuncs.com/compatible-mode/v1`**（华北2·北京） |
| region | `beijing` |
| model_id | `deepseek-v4.1-flash` |
| 来源 | `https://help.aliyun.com/zh/model-studio/base-url`（本机 2026-10-06 抓取，HTTP 200） |

⚠️ **旧命令使用的 `token-plan.cn-beijing.maas.aliyuncs.com` 已被代码明确拒绝。**
那属 **Token Plan（个人版）**，官方条款限定「仅限在编程工具和智能体工具中**交互式使用，
不可用于自动化脚本、自定义应用程序后端或任何非交互式批量调用场景**」，并警告越界使用
可能导致订阅被暂停或 Key 被封禁。本 benchmark 是批量脚本，属被禁情形。
`client.BLOCKED_ENDPOINTS` 在**发出任何请求之前**抛出 `EndpointNotPermitted`。
**不提供伪装交互或重命名脚本的绕过方式。**
详见 `../08_deliverables/channel_migration.json`。

⚠️ 该模型默认**开启思考**，会把整个 `max_tokens` 预算耗在推理上并返回**空正文**。
本项目已在 `run_config.json` 的 `generation.extra_body` 里设置
`enable_thinking: false`。**关闭思考属于推理条件变更**，开启思考的试验必须另开一张表。
详见 `../08_deliverables/benchmark_protocol.md` 第 11 节。

> 「Batch 批量推理」是另一个产品，该模型**不支持**它；这**不影响**普通
> OpenAI 兼容接口的脚本调用，二者是两回事。

### 2.1 需要你配置的四项

**provider 是显式配置项**，`config/model_config.json` 的 `active.provider` 决定用哪条通道：

| provider | 用途 | 凭据变量 | 默认 endpoint |
| --- | --- | --- | --- |
| `bailian`（默认） | 普通百炼模型 API，按量计费 | `BMB_API_KEY` → `DASHSCOPE_API_KEY` | `https://dashscope.aliyuncs.com/compatible-mode/v1` |
| `token_plan` | Token Plan 个人版，订阅套餐 | **`BMB_TOKEN_PLAN_API_KEY`** | `https://token-plan.cn-beijing.maas.aliyuncs.com/compatible-mode/v1` |

> `token_plan` 由**用户明确选择**。本项目**不声称该用途已获平台许可**：平台明示套餐
> 「仅限在编程工具和智能体工具中交互式使用，不可用于自动化脚本……非交互式批量调用」。
> 原文与来源见 `client.TOKEN_PLAN_TERMS_CN` / `channel_migration.json`，并会写入每个运行 manifest。
> 客户端**不搜索其它工具的凭据文件、不复制密钥、不读取 Claude Code 的变量**。

| 配置项 | 放在哪 | 说明 |
| --- | --- | --- |
| 对应的 Key 变量 | 环境变量 | `bailian` 用 `BMB_API_KEY`；`token_plan` 用 `BMB_TOKEN_PLAN_API_KEY`。代码不写死、不落盘、不进报告。 |
| `BMB_REGION` | 环境变量或 `config/model_config.json` | Key **按地域绑定**，必须与创建 Key 的地域一致。 |
| `BMB_BASE_URL` | 同上 | 该 provider 的 base_url，见上表。 |
| `BMB_MODEL_ID` | 同上 | 精确模型 ID，从控制台模型列表复制。 |

地域与 base_url（来自阿里云官方文档，详见 `config/model_config.json`）：

| `BMB_REGION` | 地域 | 现行 base_url | 旧域名（文档称仍可用） | 有免费额度 |
| --- | --- | --- | --- | --- |
| `beijing` | 华北2（北京） | `https://{WorkspaceId}.cn-beijing.maas.aliyuncs.com/compatible-mode/v1` | `https://dashscope.aliyuncs.com/compatible-mode/v1` | **是** |
| `singapore` | 新加坡 | `https://{WorkspaceId}.ap-southeast-1.maas.aliyuncs.com/compatible-mode/v1` | `https://dashscope-intl.aliyuncs.com/compatible-mode/v1` | 否 |
| `hongkong` | 中国香港 | `https://{WorkspaceId}.cn-hongkong.maas.aliyuncs.com/compatible-mode/v1` | `https://cn-hongkong.dashscope.aliyuncs.com/compatible-mode/v1` | 否 |
| `us_virginia` | 美国（弗吉尼亚） | `https://{WorkspaceId}.us-east-1.maas.aliyuncs.com/compatible-mode/v1` | `https://dashscope-us.aliyuncs.com/compatible-mode/v1` | 否 |
| `japan_tokyo` | 日本（东京） | `https://{WorkspaceId}.ap-northeast-1.maas.aliyuncs.com/compatible-mode/v1` | 无 | 否 |

`{WorkspaceId}` 是业务空间 ID（控制台「工作空间详情」页）。填不上时用旧域名即可。

文档另外说明：**只有北京地域有免费额度**，每个模型独立、通常 100 万 token、
有效期 90 天；快照版本与 latest 算不同模型、额度不通用。

本机当前**没有任何百炼配置**。已核实的候选模型 ID（`config/model_config.json`
`model_candidates_text`）包括 `qwen3.7-plus`、`qwen3.8-flash`、`qwen3.7-flash`、
`qwen3.8-max`、`qwen-plus`、`qwen-flash`、`qwen-turbo`。**模型目录变动频繁，
运行前请以控制台为准**，不要假设文档里的 ID 仍然存在。

### 2.2 本地配置方法（二选一）

**方式 A：只对当前终端生效（推荐，不落盘）**

```powershell
$env:BMB_API_KEY  = "<在百炼控制台为普通模型 API 创建的 Key>"
$env:BMB_REGION   = "beijing"
$env:BMB_BASE_URL = "https://dashscope.aliyuncs.com/compatible-mode/v1"
$env:BMB_MODEL_ID = "deepseek-v4.1-flash"   # 换成你控制台上的确切 ID
& $PY run_benchmark.py smoke
```

不要把 Key 贴进聊天；关闭这个终端后配置即消失。

关闭这个终端后配置即消失。

**方式 B：写成本地文件，重复运行时省事**

```powershell
Copy-Item config/local.env.example config/local.env
notepad config/local.env     # 填入四个值
& $PY run_benchmark.py smoke
```

`config/local.env` 已在 `.gitignore` 中；`run_benchmark.py` 只在启动时读入环境变量，
不会把它写进任何日志、运行记录或报告。仓库里只保留不含密钥的
`config/local.env.example`。

### 2.3 第一次真实调用建议的顺序

```powershell
# 1. 先用一种表示、一个任务、极小样本确认打通（Our Code 优先）
& $PY run_benchmark.py smoke --representations our_code --tasks t3_atom_counts --run-id probe-01

# 2. 确认原始响应、解析、评分都正常后，再放开三种表示
& $PY run_benchmark.py smoke --run-id smoke-01
```

---

## 3. 本轮任务与指标

| ID | 任务 | 模型输出 | 主指标 | 状态 |
| --- | --- | --- | --- | --- |
| T1 | Identity 同一性 | `{"same_molecule": bool}` | Accuracy（+Balanced Acc、混淆矩阵） | 已实现，待测 |
| T2 | 子结构识别（官能团） | 每标签 `true/false` | Macro-F1（+全样本惩罚版） | 已实现，待测 |
| T3 | 结构理解（重原子计数） | 每元素整数 | Exact-match Accuracy | 已实现，待测 |
| T4 | **QED 估计（本项目新增）** | `{"qed": 0..1}` | MAE（+Spearman、阈值成功率） | 已实现，待测 |
| T5 | 基础描述符（HBD） | `{"hbd": int}` | Exact-match Accuracy | 已实现，待测 |

T4 的 QED 是**本项目新增任务**，不是参考论文的任务；论文的八个任务里没有 QED 回归。
T5 的 HBA / logP / TPSA 已实现但本轮**未启用**，报告中标「待测」，不混进主表。

---

## 4. 目录说明

```
07_benchmark/
├── run_benchmark.py        唯一 CLI 入口
├── bmb/                    评测模块（不改核心类）
│   ├── paths.py            路径
│   ├── core_loader.py      载入教授原版核心，并校验它与基线逐字节相同
│   ├── reps.py             表示适配器 + 注册表
│   ├── dataset.py          ChEBI-20 → 样本池 → dev/test → manifest
│   ├── labels.py           参考标签（RDKit oracle）
│   ├── tasks.py            题目、prompt、答案解析
│   ├── scoring.py          评分规则与 cluster bootstrap
│   ├── client.py           百炼 OpenAI 兼容客户端（密钥仅从环境变量读）
│   ├── runner.py           运行循环、断点续跑、重新评分
│   └── report.py           主表、配套表、progress_matrix
├── config/
│   ├── model_config.json           端点与地域（不含密钥）
│   ├── functional_groups.json      固定 SMARTS 标签集
│   ├── functional_groups.frozen.json  prepare 后按实测支持度冻结的集合
│   ├── run_config.json             样本量、种子、预算、阈值
│   └── local.env.example           本地配置模板（真文件不提交）
├── data/raw/               ChEBI-20 原始分片（MolT5 版）
├── data/prepared/          pool / dev / test / dataset_manifest.json
└── runs/<run_id>/          raw_responses.jsonl、scored.jsonl、metrics.json、summary.csv
```

---

## 5. 数据来源与筛选规则

来源：**ChEBI-20**，取 MolT5 仓库的再发布版
（`https://github.com/blender-nlp/MolT5/tree/main/ChEBI-20_data`，BSD-3-Clause）。
化学结构本身来自 ChEBI。本轮使用 `validation` + `test` 两个分片，共 6,601 行。

筛选规则（在**任何模型调用之前**执行，逐条记录在 `dataset_manifest.json`）：

| 规则 | 内容 | 剔除数 |
| --- | --- | --- |
| F1 | RDKit 能解析 | 0 |
| F2 | 源 SMILES 不含立体化学标记 | 4,498 |
| F3 | Our Code 能解析 | 0 |
| F4 | Our Code 写出后可重新解析 | 0 |
| F5 | RDKit 认可 Our Code 的输出是同一分子 | 0 |
| F6 | RDKit 规范 SMILES 也能被 Our Code 读取 | 0 |
| F7 | 规范 SMILES 重新解析后不变（引用唯一） | 3 |
| F8 | SELFIES 能编码 | 1 |
| — | 同一分子重复 | 2 |

**F2 是最大的一道筛**：ChEBI-20 有 68% 的分子带立体化学。教授原型明确声明
不支持立体化学并且**拒绝**而不是丢弃，所以这些分子被**排除**，不是剥掉立体信息后
继续用——那样会改变化合物身份，等于先丢信息再声称还是原样本。

结果：**2,097** 个三种表示都能无损表达的共同分子，dev 524 / test 1,573，
分子级划分（同一分子的所有变体、Identity 题对两端都在同一 split）。

---

## 6. 公平性规则（已固化在代码里）

- 三种表示用**同一份分子清单、同一题 ID**，逐题配对比较。
- prompt 除「表示说明一句话」与载荷外完全相同；每个 prompt 存 hash。
- 模型全程**不调用工具、不执行 Our Code、不调用 RDKit**。
- 每题**独立会话**，不共享其它表示的回答；调用顺序用固定种子打乱。
- **重试只针对传输故障**（超时、连接错误、429/5xx）。看到错答案后重试样一次是禁止的。
- 无效应答在全样本指标里记失败；**不裁剪、不当作误差 0、不默认为某个正确类别**。
- 未跑完的格子标「待测」，不排名，不因为样本少就宣称某表示更优。
- 置信区间用 **cluster bootstrap**（按母分子重采样），F1 与相关系数在每次重采样内重算。

---

## 7. 许可证与引用注意

- `arunraja-hub/molrepbench` 仓库（论文作者公开代码）在核查时**未声明任何许可证**
  （截至 commit `21d36941cddf`）。本项目**未复制其代码或数据**，只借鉴实验设计
  （表示 × 任务 × 模型、同题配对、客观标签、按任务选指标），并按论文本身引用。
- 数据用 MolT5 的 BSD-3-Clause 再发布版，不用 MolRepBench 的
  `benchmark_data_small.zip`。
- 论文正文对 logP/TPSA 用 Spearman，附录另一处表述不同；本项目**事先声明**
  回归任务用 MAE + Spearman，不复制含混表述。
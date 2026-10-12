# 教授原始基准（只读，不得修改）

本目录是本轮研究的**唯一基准来源**。所有复现都针对这里的未修改副本进行。

## 1. 文件清单与哈希（SHA-256）

| 文件 | SHA-256 | 来源 | 说明 |
| --- | --- | --- | --- |
| `molecule_standalone.py` | `55d314aa312707b498536d339ff59dc89bbbe4956cab54a74a4746cd385af57b` | 教授提供，复制自项目根目录同名文件 | 教授原型完整源码（490 行） |
| `貼上的 Markdown (1).md` | 见 `BASELINE_HASHES.txt` | 教授提供，复制自项目根目录 | 教授聊天内容（研究说明基准） |
| `raw.txt` | 见 `BASELINE_HASHES.txt` | 教授提供 | **历史**测试日志（16 tests OK），非本轮结果 |
| `verification_results.json` | 见 `BASELINE_HASHES.txt` | 教授提供 | **历史** RDKit 对照报告（4991 molecules），非本轮结果 |

## 2. 与教授旧报告的哈希核对

`verification_results.json` 内字段 `module_sha256` = `55d314aa312707b498536d339ff59dc89bbbe4956cab54a74a4746cd385af57b`。

本轮实测 `molecule_standalone.py` 的 SHA-256 与该值**逐字符相同**。

**结论**：教授提供的原型文件，正是旧对照报告所验证过的那一版文件，文件本身未被改动过。

## 3. 本轮不可复现的原始实验（材料缺失）

| 教授自述的实验 | 依据文本 | 缺失材料 | 本轮处理 |
| --- | --- | --- | --- |
| 16 个测试方法全部通过 | `貼上的 Markdown (1).md` 第 177 行；`raw.txt` | `test_standalone.py` **不在工作区** | 无法运行原套件；本轮自行重写测试并单独报告，**不复用 raw.txt 的结论** |
| 60 个具名 SMILES 用例、200 组原子/键序排列 | 同上 | 同上 | 同上 |
| 4999 条 NCI 语料、4991 条与 RDKit canonical SMILES 一致 | `貼上的 Markdown (1).md` 第 179 行；`verification_results.json` | `first_5K.smi` **不在工作区**；且本轮禁止调用 RDKit | **当前材料下未完成原始完整实验复现**。本轮不引用该数字作为成果 |

## 4. 保护措施

- 活动项目根目录的 4 个教授附件保持**原位、原样、未修改**。
- 归档目录 `c:/Users/Asus/Desktop/Design Molecule Class_archive_20260928_175738/protected_materials/` 另存一份副本。
- 本轮所有改进只发生在 `research/01_implementation/`，不写回本目录。

## 5. 注意：本项目内存在与本任务无关的旧环境

项目根目录下的 `.venv310/` 是上一轮为安装 `rdkit==2026.3.6` 建立的虚拟环境（见归档项目的 `requirements.txt`）。

本轮**不使用、不激活**该环境，所有复现与测试使用系统解释器 `C:\Users\Asus\AppData\Local\Programs\Python\Python39\python.exe`（Python 3.9.0）。按用户要求未扩大清理范围，故该目录原地保留。

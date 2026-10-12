# 01_implementation — 交付核心（教授原版）

## 这里是什么

`molecule_core.py` 是**本轮交付的核心实现**，内容与教授原版
`../00_baseline/molecule_standalone.py` **逐字节相同**：

```
SHA-256  55d314aa312707b498536d339ff59dc89bbbe4956cab54a74a4746cd385af57b
```

这个值同时等于项目根目录 `molecule_standalone.py` 的哈希，也等于教授自带
`verification_results.json` 里记录的 `module_sha256`。

保留两份相同文件的原因：`00_baseline/` 是**不可修改的历史基准记录**，
本目录是**对外交付与 import 的工作路径**。两者相同这件事由测试强制保证，
不靠人工记着。

## 本轮撤销了什么

上一轮在核心类上做过四处自行扩展，都不属于教授原版范围，本轮已全部撤销：

| 撤销项 | 位置 |
| --- | --- |
| `Molecule.is_isomorphic()` 回溯同构判定 | `Molecule` |
| `label_of()` 与 `ATOM_LABEL_FIELDS`（仅为 `is_isomorphic` 服务） | 模块级 |
| 环闭合「可回收标签池」`_closure_text()` / `_next_label()` | `write_smiles` |
| `to_code()` 的长度自检与 `ast.parse` 自检 | `Molecule.to_code` |

撤销前的完整内容保存在 `../99_archive_rollback_20261005/`。

**为什么连 `is_isomorphic` 也要撤**：Identity 任务的正确答案必须由**外部固定化学工具**
（本项目用 RDKit）产生。如果核心类自己带同构判定，再拿它给模型打分，等于让被测对象
自己出标准答案。判定真值的代码现在在 `../07_benchmark/bmb/labels.py`，
不在这三个类里。

## 保留了什么

教授原有的全部必要功能都在：SMILES 解析、写出、图记录校验、邻居/键/最短路查询、
原子与键编辑、总电荷、构造代码读写，以及 `from_code()` 的资源限制
（`MAX_TEXT = 100000`、AST 节点上限 40000、`MAX_ATOMS = 10000`）。
「只保留三个基本类」指的是保持原版的**实现范围**，不是把解析器等支撑代码删掉。

## 谁在用

- 核心测试：`../02_tests/test_baseline_core.py`（55 项，含 2 项防回归）
- benchmark：`../07_benchmark/bmb/core_loader.py` 按路径加载它，
  并在加载前校验它与基线逐字节相同

## 运行

核心本身**只用 Python 标准库**，用系统解释器即可（不需要项目 venv）：

```powershell
cd "../02_tests"
python -S -B -m unittest discover -s . -p "test_*.py"
```

benchmark 模块需要 rdkit / selfies，必须用 `.venv` 的解释器，见
`../07_benchmark/README.md`。
# FICPF GitHub 精简发布版

FICPF（Frozen Inductively Calibrated Pathomic Fusion）是一种面向病理-组学生存分析的冻结式归纳校准融合方法。方法保留已经训练好的 VLSA-CONCH 病理预测器，只在每个 outer training fold 内拟合低容量 RNA-seq/CNV PCA-Cox 分支；随后将两种预测映射到训练折经验分位数空间，并按固定的病理/组学权重 `0.80/0.20` 融合。

![FICPF 方法图](figures/ficpf_20260821/ficpf_method_overview.png)

## 锁定结果

- matched VLSA 五队列宏平均 C-index：`0.695782`。
- FICPF 五队列宏平均 C-index：`0.705017`。
- 配对宏平均增益：`+0.009235`，95% bootstrap CI 为 `[0.0016, 0.0169]`。
- 数据集：BLCA、BRCA、GBMLGG、LUAD、UCEC，均采用患者级 outer 5-fold。
- VLSA 病理主干使用 `42、3404、3405` 三个完整训练种子进行稳定性分析，15 个“种子 × 队列”点估计全部正向。

正式主结果仍是锁定的 seed 42 结果，多种子实验只用于稳定性验证，不反向替换主表。

## 为什么小于 25 MB

该目录是面向 GitHub 的白名单精简包。它保留：

- FICPF 全部分析和制图代码；
- 固定数据划分、生存表和逐患者锁定预测；
- 主结果、消融、敏感性、鲁棒性和多种子统计结果；
- 中英文论文草稿、参考文献和正式图件；
- 环境信息与逐文件 SHA256 清单。

它不包含原始 WSI、CONCH 权重、patch 特征、VLSA checkpoint、大型组学表、W&B、日志、缓存和失败实验。大型组学表的文件名、放置位置和 SHA256 见 [DATA.md](DATA.md)。

## 安装

```bash
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows: .venv\Scripts\activate
python -m pip install -r requirements.txt
```

## 无大数据校验

```bash
python tools/verify_ficpf_release.py
```

该命令不需要被省略的组学表，会检查文件清单、包体积、主结果、多种子结果和禁止混入的训练产物。

## 复算边界

1. 使用 `result/` 中锁定逐患者预测进行汇总复核和制图，可直接运行。
2. 完整组件消融、外部基线统一比较、重新拟合 PCA-Cox 或组学损坏实验，需要按 [DATA.md](DATA.md) 补齐组学表。
3. 从原始 WSI 重训 VLSA 还需要上游 VLSA 代码、CONCH 权重和 patch 特征，不属于本精简分析仓库。

完整英文说明见 [README.md](README.md)，许可证和第三方边界见 [NOTICE.md](NOTICE.md)。

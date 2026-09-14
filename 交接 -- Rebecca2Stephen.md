# 交接 -- Rebecca2Stephen

> 交接人：Rebecca  
> 接收人：Stephen  
> 日期：2026-09-14  
> 项目方向：AAV capsid 变体的包装效率与组织嗜性预测

## 1. 背景与参考仓库

本次方案参考以下两个仓库：

- [Fit4Function](https://github.com/vector-engineering/fit4function)：AAV capsid 多功能筛选与预测。公开环境基于 Python 3.9、TensorFlow 2.13、Keras 2.12，并提供 `Final_model_24K_v1.h5` 预训练模型。数据说明中包含约 5 万、10 万、24 万规模的变体数据集。
- [LIFS-Comp_AAV9](https://github.com/stephenovo/LIFS-Comp_AAV9)：面向 AAV9/SMA 的计算筛选框架，已经采用 Python 3.11、7-mer one-hot、Ridge、Random Forest，并预留 LightGBM 和 PyTorch 多任务模型。

## 2. 推荐结论

不建议直接复刻 Fit4Function 的旧 TensorFlow 环境作为新项目主线。建议采用：

```text
Python 3.11
PyTorch 2.2+
ESM-2（先冻结参数提取 embedding）
LightGBM + scikit-learn
Polars/Pandas + PyArrow + Parquet
MLflow 或 Weights & Biases
Git + DVC + Docker
```

核心策略是：先建立可解释的传统机器学习基线，再验证蛋白语言模型 embedding 是否在严格的序列隔离测试集上带来真实提升。

## 3. 编码器和模型路线

### 阶段 A：可解释基线

1. 7-mer amino-acid one-hot，输入维度约为 `7 × 20 = 140`。
2. 增加氨基酸理化性质：电荷、疏水性、体积、芳香性等。
3. 分别训练 Ridge、Random Forest、LightGBM。

这一步不需要 GPU，可以快速检查标签质量、数据泄漏和任务可预测性。

### 阶段 B：冻结 ESM-2 embedding

推荐顺序：

```text
ESM-2 8M：快速原型
ESM-2 35M：正式起点，性能/成本平衡较好
ESM-2 150M：数据规模明显扩大后再考虑
```

不要只把 7 个插入氨基酸单独输入蛋白模型。优先使用：

- 插入位点上下游 20–40 aa 窗口；或
- 完整 AAV9 VP1 序列；
- 同时保留 7-mer one-hot 作为显式特征。

ESM embedding 后接 LightGBM 或小型 MLP，通常比直接微调大模型更适合目前的数据规模。

### 阶段 C：多任务预测

可对以下 endpoint 建立独立输出头：

- packaging/production
- brain
- spinal cord
- liver
- heart/kidney 等 off-target

每个任务独立处理缺失标签。包装效率建议作为硬门槛或单独任务；综合评分和 Pareto 排名放在预测之后计算，不要一开始把最终评分直接作为训练标签。

## 4. 数据和验证要求

建议统一为 canonical schema：

```text
variant_id
aa_sequence
insertion_sequence
library_id
batch_id
animal_id
replicate_id
packaging_score
brain_score
spinal_cord_score
liver_score
heart_score
kidney_score
missing_label_flags
```

重要要求：

1. 不要随机切分相似序列。优先按 library、batch、动物或序列家族做 sequence-aware split。
2. 训练、验证、测试集之间要避免同一变体或高度相似变体泄漏。
3. 报告 MAE、RMSE、Spearman/Pearson、命中率和置信区间，而不只报告单一 R²。
4. 最终候选要同时检查包装门槛、预测不确定性、序列距离、变体多样性和结构可行性。

## 5. 推荐服务器

### 推荐起步配置

```text
Ubuntu 22.04/24.04
16 vCPU
64 GB RAM
NVIDIA L4 24 GB 或 RTX 4090 24 GB
2 TB NVMe
```

这套配置足以完成基线训练、ESM-2 embedding 提取、LightGBM 和小型 PyTorch 多任务模型。

### 其他配置

| 用途 | 配置 |
|---|---|
| 仅基线和数据清洗 | 8–16 vCPU，32–64 GB RAM，500 GB–1 TB NVMe，无 GPU |
| 常规正式训练 | 16 vCPU，64 GB RAM，L4/4090 24 GB，1–2 TB NVMe |
| 大规模 embedding 或模型微调 | 24–32 vCPU，128 GB RAM，A100 40/80 GB，2–4 TB NVMe |

目前不建议一开始购买多张 A100。先确认数据覆盖率和严格测试集上的增益，再决定是否扩容。

## 6. 容量估算

处理后的表格通常不是主要存储压力：

| 变体数量 | 处理后表格 | ESM embedding 粗略规模 |
|---:|---:|---:|
| 10 万 | 0.1–1 GB | 0.2–1 GB |
| 100 万 | 1–10 GB | 2–10 GB |
| 500 万 | 10–50 GB | 10–50 GB |

实际规划：

- 1 TB：最低可用；
- 2 TB：推荐，适合训练和保留中间结果；
- 4–8 TB：需要保留原始 FASTQ/BAM 和多批次测序数据时使用。

建议 embedding 保存为 float16，数据集保存为 Parquet，并把原始数据、中间数据、处理数据和模型 artifact 分开管理。

## 7. 建议执行顺序

1. 完成原始数据审计和 canonical schema 映射。
2. 用 7-mer one-hot 建立 Ridge/RF/LightGBM 基线。
3. 做 sequence-aware train/validation/test split。
4. 提取 ESM-2 35M embedding，与基线进行同一测试集比较。
5. 训练多任务 MLP，并加入 ensemble 或不确定性估计。
6. 用包装门槛、CNS reward、liver/off-target penalty 做候选排序和 Pareto 分组。
7. 形成候选短名单后，再进入体外和动物实验验证。

## 8. 交接边界

这些模型输出的是 AAV 包装和组织分布的计算代理指标，不等同于人体疗效、安全性或临床可行性。尤其是“降低预测肝脏富集”不能直接解释为“降低肝毒性”，最终结论必须依赖实验验证。

## 9. 待 Stephen 确认

- 最终数据是否包含原始 FASTQ/BAM，还是只有整理后的变体标签表；
- 目标是复现 Fit4Function，还是优先完成 AAV9/SMA 的新候选筛选；
- 可接受的服务器来源和预算；
- 是否需要把训练流程包装成可重复运行的 CLI/Docker pipeline。


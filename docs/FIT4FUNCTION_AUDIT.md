# Fit4Function 深度数据可行性审计

> 审计日期：2026-09-14；原始数据验证更新：2026-09-15
> 项目：LIFS-Comp_AAV9  
> 结论等级：**黄灯转绿 - 科学路线可行，但公开处理表不能直接训练 CNS 多任务模型，必须先重建 SRA 原始测序数据。**

## 1. 最先看这一页：我们到底能不能做？

可以做，但不是“下载一个 CSV 就开始训练”。

Fit4Function 官方 GitHub 的处理后数据分成两类：

1. **带 7-mer 序列的数据**：包装、鼠肝、HepG2、THLE-2 等标签可以立即训练；
2. **多器官数据**：脑、脊髓、肝、心、肾等表只有匿名 `SequenceID`，没有 7-mer，也没有公开映射表，因此不能直接训练序列模型。

真正让项目重新可行的是 NCBI SRA。BioProject `PRJNA1131359` 已公开 270 个实验、约 120.32 GB 测序数据，包含 Fit4Function 库在四只小鼠的脑、脊髓、肝、心、肾等原始 reads。[^1][^5] 我们先用一条脊髓技术重复验证 R1 → 21 nt → 7-mer，再完整重建 12 个 Liver 与 6 个 prod2/prod3 Virus DNA runs。以 prod2 为分母得到的 100K 肝标签与官方 `Liver` 标签达到 Pearson `r = 0.9778`，因此原始数据重建路线已经得到外部标签验收。

### 最终判断

| 模块 | 只用 GitHub 处理表 | 加上 SRA 原始 reads | 当前判断 |
| --- | --- | --- | --- |
| `F_pack` 包装 | 可训练 | 可复核 | 绿灯 |
| `F_liv` 小鼠肝 | 可训练，100K 条带序列标签 | 可重建并交叉验证 | 绿灯 |
| 人肝细胞辅助惩罚 | 可训练 | 可复核 | 绿灯 |
| `F_CNS` 脑/脊髓 | **不可直接训练** | 可从 raw reads 重建 | 黄灯转绿 |
| `F_off` 心/肾 | **不可直接训练** | 可从 raw reads 重建 | 黄灯转绿 |
| 人运动神经元特异性 | 没有数据 | SRA 也没有 | 红灯，不能声称 |
| 肝毒性降低 | 没有毒性终点 | SRA 也没有 | 红灯，只能说“肝分布代理下降” |

## 2. Fit4Function 实际做了什么

论文研究的是 AAV9 衣壳表面 588/589 位置之间插入 7 个氨基酸的变体。7-mer 理论空间为 `20^7 = 12.8 亿`，不可能全部做实验，因此作者先学习“能不能生产/包装”，再从预测可生产空间均匀抽取 24 万条，建立 Fit4Function 库，最后分别训练不同功能模型。[^1]

需要特别注意：论文实验骨架写作 **AAV9 (K449R)**。因此这些数据不是对任意 AAV9 骨架都自动成立，更不能直接等同于 Zolgensma 的临床衣壳背景。我们的候选定义必须把“Fit4Function AAV9(K449R) 588/589 7-mer 插入背景”写清楚。

### 关键实验设计

| 项目 | 论文设置 | 对我们的影响 |
| --- | --- | --- |
| 设计单元 | 7-mer 氨基酸序列 | 当前 one-hot / LSTM 输入设定正确 |
| 插入位置 | VP1 588 与 589 之间 | 必须固定，不外推到其他位点 |
| 包装建模库 | 两个约 74.5K 库，约 10K 共享 | 可做独立库外测试 |
| Fit4Function 库 | 约 240K 预测可生产变体 | 多器官筛选的母库 |
| 小鼠 | 4 只成年雌性 C57BL/6J | 样本量与性别/年龄外推有限 |
| 给药 | 静脉、每鼠 `1×10^12 vg` | 是全身分布筛选，不是 SMA 疗效实验 |
| 取材 | 给药后 2 小时 | 测的是早期载体基因组分布 |
| 器官 | 血清、肝、脾、肾、心、肺、脊髓、脑 | 可构建我们的 CNS/肝/脱靶代理 |
| 功能分数 | `log2(器官 RPM / 病毒库 RPM)` | 需要病毒库作为分母重建 |

论文使用 150 bp 单端 R1；将 reads 比对到带 21 个 `N` 的 AAV9 短参考序列，提取 21 nt 插入，并过滤插入区任一碱基 QScore < 20、含 `N`、有 indel 或不属于合成库的 reads。之后按样本计算 RPM、平均技术重复，再相对生产后病毒库计算 log2 enrichment。[^1]

## 3. 我们审计了哪些官方来源

| 来源 | 版本/规模 | 用途 | 结果 |
| --- | --- | --- | --- |
| 论文正文 | Nature Communications 2024 | 方法、实验边界、模型设置 | 已逐段核对 |
| 官方 GitHub | 最新提交 `6bfc2ebfe4` | 处理表、notebook、生产模型 | 已完整克隆和哈希 |
| Zenodo 软件归档 | v1.0.0，73.9 MB | 冻结版代码/数据 | 与 GitHub v1.0.0 文件清单一致[^3] |
| Zenodo Source Data | 175.0 MB | 复现论文图 | 是画图源数据，不补齐匿名序列映射[^4] |
| NCBI SRA | 270 runs，120.32 GB | 原始 NGS reads | 多器官重建的关键来源[^5] |
| 补充材料 | 10 页 | 重复性与补充验证 | 已下载并核对[^6] |

官方代码采用 BSD 3-Clause 许可证，可以在保留版权和许可文本的前提下修改和再分发。[^8]

## 4. GitHub 处理数据：逐文件审计

所有机器可读结果见：

- [`audit_data/fit4function_release_audit.json`](audit_data/fit4function_release_audit.json)
- [`audit_data/fit4function_sra_manifest.csv`](audit_data/fit4function_sra_manifest.csv)
- [`audit_data/fit4function_sra_target_manifest.csv`](audit_data/fit4function_sra_target_manifest.csv)
- [`audit_data/fit4function_sra_summary.json`](audit_data/fit4function_sra_summary.json)

### 4.1 包装数据

| 文件 | 行数 | 唯一合法 7-mer | 备注 |
| --- | ---: | ---: | --- |
| `modeling_library_production_fitness.csv` | 74,464 | 74,464 | 两个 codon replicate + 聚合生产值 |
| `assessment_library_production_fitness.csv` | 74,467 | 74,467 | 可作独立库外测试 |

两个文件的实际氨基酸序列交集为 9,993 条，接近论文所述 10K 共享集。各文件含 9,987 条标为 `Calibration` 的记录；额外交集可能来自随机序列碰撞。这个细节不影响使用，但切分时必须按**实际序列交集**去重，不能只相信标签。

我们对 modeling library 的有效 `Designed` 生产值取 log2，并拟合两个高斯分量：

| 分量 | 均值 | 标准差 | 权重 |
| --- | ---: | ---: | ---: |
| 低生产/非适应分量 | -7.323 | 2.068 | 0.559 |
| 可生产分量 | 0.233 | 1.897 | 0.441 |

两分量后验概率相等的审计边界约为 **log2 enrichment = -3.304**，对应原始比例约 0.101。有效记录中约 43.7% 高于该边界，约 26K 条，与论文图中的 production-fit 数量相符。这个数是可复现的初始估计，不应未经验证直接写死为最终硬门槛。

### 4.2 带序列的 100K 功能筛选表

`fit4function_library_screens.csv` 有 100,000 条互不重复且格式合法的 7-mer，包含：

- `Production1`, `Production2`
- `Liver`
- `HepG2_bind`, `HepG2_tr`
- `THLE_bind`, `THLE_tr`

因此包装、鼠肝和人肝细胞代理可以立即做序列建模。缺失/无穷值不是小问题，尤其转导列：

| 标签 | 有限值 | 可用率 | 主要问题 |
| --- | ---: | ---: | --- |
| Production2 | 98,568 | 98.6% | 少量缺失 |
| Liver | 98,621 | 98.6% | 少量缺失和 ±Inf |
| HepG2 binding | 76,867 | 76.9% | 大量未检出 |
| HepG2 transduction | 32,924 | 32.9% | 约 65K 为 `-Inf` |
| THLE binding | 91,509 | 91.5% | 缺失较少 |
| THLE transduction | 54,096 | 54.1% | 约 38.7K 为 `-Inf` |

不能简单把所有 `-Inf` 删除后宣称模型覆盖完整分布。至少要比较三种处理：未检出作为左删失、设定检测下限、以及只对有限值回归，并在报告中说明差异。

### 4.3 多器官测量表的致命缺口

`fit4function_library_invivo.xlsx` 有 8 个 sheet，每个 50,000 行、四只动物。但其 ID 形如：

```text
Sequence_Brain_1
Sequence_SpinalCord_1
Sequence_Liver_1
```

工作簿没有 `AA` 列、没有隐藏 sheet、没有 defined name，也没有 ID-to-AA lookup。更重要的是，将数字后缀强行对齐后，各器官均值的相关系数几乎都在 0 附近，说明不同 sheet 至少经过独立排列，不能把 `_1` 当作同一条 7-mer。

补充材料的 Supplementary Fig. 4 明确展示了同一变体在不同器官之间存在成块的正负相关结构；这与强行按匿名数字后缀拼接得到的近零相关明显不符，进一步证明公开 workbook 的跨 sheet 行号不是共享序列键。[^6]

所以，**不能把 Brain 第 1 行与 Liver 第 1 行拼起来，也不能拿这些匿名表训练序列模型。** 这是本次审计最重要的发现。

### 4.4 多器官数据质量

匿名表仍可用于检验论文的重复性和预测性能：

| 器官 | 四只动物两两 Pearson r 中位数 | 范围 | 官方测量-预测 r |
| --- | ---: | ---: | ---: |
| Serum | 0.771 | 0.723–0.802 | 0.690 |
| Liver | 0.843 | 0.788–0.887 | 0.822 |
| Spleen | 0.817 | 0.790–0.842 | 0.744 |
| Kidney | 0.765 | 0.725–0.793 | 0.707 |
| Heart | 0.526 | 0.442–0.613 | 0.619 |
| Lung | 0.474 | 0.222–0.752 | 0.575 |
| Spinal cord | 0.661 | 0.625–0.739 | 0.646 |
| Brain | 0.589 | 0.551–0.647 | 0.612 |

这正好解释为什么肝模型比脑/脊髓模型更稳定：模型上限受实验重复性限制。对我们的 SMA 故事，`F_CNS` 是可建模的，但噪声会明显高于 `F_liv`，最终必须报告动物留出不确定性。

## 5. SRA 原始数据审计

NCBI 当前公开 270 个 runs、41 个 BioSamples、约 318 Gbases/0.12 TB。[^5] 提交者使用的库代号可对应为：

| SRA library | 推断用途 | runs | SRA 记录字节数 |
| --- | --- | ---: | ---: |
| Hammerhead | Fit4Function 库 | 127 | 61.49 GB |
| Multifunction | MultiFunction 验证库 | 83 | 28.34 GB |
| shark2 / shark3 | 包装 modeling / assessment 库 | 36 | 21.22 GB |
| NNK | 随机 NNK 对照 | 18 | 3.17 GB |
| Winghead | 猕猴相关库 | 6 | 6.10 GB |

`Hammerhead` 中存在完整的四只动物与技术重复：脑、脊髓、肝、脾、心、肺、肾各 `4 × 3 = 12` runs；血清每只动物 1 个 run；另有 9 个 Virus DNA runs 和 3 个 Plasmid DNA runs。

这里还发现一处元数据矛盾：官方 GitHub README 把 in-vivo workbook 描述为每个器官、每只动物平均“四个技术重复”，但论文 Methods 写的是 triplicate，SRA 文件结构也明确是多数器官每只动物 3 个技术重复。[^1][^2] 后续以论文 Methods 与原始 run manifest 为准，并把四只动物与三个技术重复分开记录。

### 只做本项目需要下载多少？

| 范围 | runs | SRA 记录规模 |
| --- | ---: | ---: |
| Brain + spinal cord | 24 | 5.96 GB |
| Liver | 12 | 5.46 GB |
| Heart + kidney | 24 | 10.44 GB |
| Virus DNA 分母 | 9 | 14.29 GB |
| **本项目核心合计** | **69** | **36.16 GB** |

无需下载整个 120 GB。若边下载边流式计数，不解压落盘，建议预留 60–100 GB 临时空间；网络会比计算更可能成为瓶颈。

## 6. 我们已经真的跑过一条脊髓 FASTQ

代表 run：`SRR29692586`，`hammerhead_biod_SpinalCord_a1_r3`。[^7]

文件通过 ENA 镜像取得，压缩后 92,386,737 bytes，MD5 与官方镜像记录一致。我们使用论文紧邻插入区的左右锚点、21 nt 长度、插入区 Q20 和标准密码子翻译做了一个**保守 pilot**：

| 指标 | 结果 |
| --- | ---: |
| 检查 reads | 1,787,085 |
| Q20 且锚点通过、无终止密码子的 reads | 1,023,545 |
| 保守规则有效率 | 57.27% |
| 有效 21-mer 种类 | 246,131 |
| 有效 7-mer 种类 | 175,699 |
| 落入公开 100K 7-mer 清单的 reads | 412,827 |
| 单个技术重复检测到的 100K 清单序列 | 70,428（70.43%） |

完整机器结果见 [`audit_data/SRR29692586_spinal_cord_pilot.json`](audit_data/SRR29692586_spinal_cord_pilot.json)。

随后我们用论文公开参数补跑了 Bowtie2 流式比对。它对同一批 1,787,085 reads 的总体
比对率为 99.77%，Q20、无插入区 indel 且无终止密码子的有效率为 59.00%；其中
419,561 reads 落入公开 100K 清单，检测到 70,742 条白名单序列。与保守锚点法相比，
白名单原始计数的 Pearson 相关系数为 0.9985，`log1p(count)` 相关系数为 0.9977。
这说明轻量 pilot 没有扭曲主要丰度结构，同时 Bowtie2 找回了锚点含少量错配的 reads。
机器结果见
[`audit_data/SRR29692586_bowtie2_pilot.json`](audit_data/SRR29692586_bowtie2_pilot.json)。

这一步证明“raw read → Bowtie2 → 21 nt → 7-mer → 公开 100K 序列”的链路可走通。
但我们仍没有完整 240K 合成库白名单，因此 59.00% 只能当作公开子集重建的管线指标，
不能写成对作者原始过滤率的精确复现。

### 6.2 Liver 外部标签验收

第一轮完整批次包含 12 个 Liver runs（4 只动物 × 3 个技术重复）以及 prod2、prod3
各 3 个 Virus DNA runs。18 个 ENA FASTQ 共 6.258 GB；下载器逐文件验证压缩字节数和
ENA MD5。Bowtie2 共处理 132,464,063 reads，其中 Liver 组 92,235,271 reads、
prod2 病毒库 23,808,806 reads、prod3 病毒库 16,419,986 reads。

我们没有先指定一个“最有利”的分母，而是同时比较：

- 全部有效 7-mer reads 或公开 100K 白名单 reads 作为 RPM 分母；
- Liver Animal 1–4、Animal 1–3 和单动物；
- prod2、prod3 或两个可用生产轮次均值作为病毒库分母。

最佳且与官方数据定义一致的组合是：白名单 RPM、四只动物均值、prod2 病毒库分母。
在 97,873 条双方均为有限值的序列上：

| 验收指标 | 结果 |
| --- | ---: |
| Pearson r | **0.9778** |
| Spearman r | **0.9822** |
| 线性校准斜率 | 0.9546 |
| 线性校准截距 | -0.0258 |
| 校准后 RMSE | 0.3092 |

“全部有效 reads”分母得到的 Pearson r 也为 0.9778，说明公开 100K 子集并没有通过
分母选择制造相关性。Animal 1–3 均值对官方标签仍有 `r = 0.9605`，但四只动物均值更高，
符合官方 `Liver` 列使用全部动物聚合的解释。prod3 单独作分母只有 `r = 0.8658`，
因此后续多器官重建固定使用 **prod2**，不需要额外下载体积约 8.78 GB 的 prod1。

完整组合结果见
[`audit_data/fit4function_liver_validation_metrics.csv`](audit_data/fit4function_liver_validation_metrics.csv)，
逐 run 质量控制见
[`audit_data/fit4function_liver_run_qc.csv`](audit_data/fit4function_liver_run_qc.csv)。

## 7. 已跑的序列基线

我们先对带序列的 100K 表进行了固定随机种子、20% 随机留出的 sanity check。它只回答“信号能否从序列中学到”，不作为最终泛化结论。

| 标签 | Ridge Pearson r | Ridge R² | Random Forest Pearson r |
| --- | ---: | ---: | ---: |
| Production2 | 0.657 | 0.432 | 0.643 |
| Liver | 0.742 | 0.550 | **0.767** |
| HepG2 binding | 0.604 | 0.365 | 未在本轮重训 |
| HepG2 transduction | 0.564 | 0.318 | 未在本轮重训 |
| THLE binding | 0.523 | 0.273 | 未在本轮重训 |
| THLE transduction | 0.589 | 0.347 | 未在本轮重训 |

结果见 [`audit_data/fit4function_baseline_metrics.csv`](audit_data/fit4function_baseline_metrics.csv) 和 [`audit_data/fit4function_random_forest_metrics.csv`](audit_data/fit4function_random_forest_metrics.csv)。简单 Ridge 已有明显信号，随机森林对 Liver 再提升，说明开始做模型比较是合理的；但最终必须改用序列距离/聚类切分和动物留出。

Supplementary Fig. 9 还显示，一阶残基效应不能解释全部预测；加入二阶相互作用后，对原模型预测的解释度明显提高。[^6] 因此 Ridge 适合作为可解释下限，但最终模型至少需要显式成对特征、树模型或能学习位置间相互作用的 MLP/LSTM。

包装模型还做了更严格的独立库外 sanity check：固定抽取 modeling library 的 24,000 条训练，排除与其重叠的序列后，在 assessment library 的 57,724 条序列上测试。Ridge 得到 `r = 0.767`，Random Forest 得到 `r = 0.780`。结果见 [`audit_data/fit4function_production_generalization.csv`](audit_data/fit4function_production_generalization.csv)。这说明我们的轻量实现有泛化信号，但尚未追平论文专门为该任务训练的 LSTM；当前应把它当作管线基准，而不是最终包装头。

## 8. 对 SMA 项目最重要的科学边界

### 8.1 这里的“脑/脊髓”是什么

它是给药后 2 小时从器官 DNA 中测得的载体基因组分布。它不是：

- 人运动神经元转导；
- 脊髓前角运动神经元特异性；
- SMN1 表达；
- 运动功能恢复；
- 长期疗效。

因此项目名称和摘要应该使用“**小鼠 CNS 早期生物分布代理**”，不能直接写“运动神经元靶向模型”。

### 8.2 `F_CNS` 不宜一开始就把脑和脊髓混成一个标签

建议先分别训练：

```text
F_brain(sequence)
F_spinal(sequence)
```

通过动物留出测试后，再将标准化预测等权合成为展示用 `F_CNS`。这样可以发现“只进脑但不进脊髓”的候选，也不会让一个高噪声器官掩盖另一个器官。

### 8.3 肝脏下降不等于肝毒性下降

`F_liv` 是相对病毒库的肝 DNA 富集。它可用于表达“预测肝分布代理较低”，但不能直接推出肝酶下降、肝衰竭风险下降或临床剂量下降。

### 8.4 AAV9(K449R) 背景不能省略

模型输入虽然只有 7-mer，但实验标签来自固定骨架。任何候选都应写成“在该实验骨架背景下的计算候选”，不能把插入肽当作脱离衣壳环境的通用靶向模块。

## 9. 建议采用的重建与训练管线

```text
SRA manifest
    ↓ 只选 Hammerhead：5 个目标器官 + Virus DNA
FASTQ 校验与流式读取
    ↓
Bowtie2 比对论文短参考序列
    ↓
提取 21 nt，Q20 / N / indel / stop 过滤
    ↓
翻译 7-mer，只保留公开 100K 清单
    ↓
每 run 的 raw counts → RPM
    ↓
技术重复平均 → 每只动物、每个器官
    ↓
log2(organ RPM / virus RPM)
    ↓
用公开 Liver 标签做重建一致性检验
    ↓
Animal 1–3 训练，Animal 4 严格盲测
    ↓
F_brain / F_spinal / F_liv / F_heart / F_kidney
    ↓
包装硬门槛 → 帕累托、多样性、不确定性
```

### 病毒库分母尚需验证

SRA 中有 3 组 production、每组 3 个技术重复，共 9 个 Virus DNA runs。公开材料没有在文件名层面明确指出每只动物对应哪一组生产批次。不能拍脑袋选分母。

我们会透明比较以下候选：

1. 各 production 组三重复平均；
2. 九个 runs 总体平均；
3. 若元数据或实验批次能建立匹配，则按批次配对。

使用公开 100K `Liver` 列作为外部校准：选择能复现公开肝标签且在重复性上合理的规则，并把选择过程记录下来。这个校准只确定数据处理，不用于调高模型测试分数。

## 10. Go / No-Go 验收条件

在开始大规模模型训练前，数据重建至少应满足：

| 检查 | 建议门槛 | 失败意味着什么 |
| --- | --- | --- |
| FASTQ 完整性 | MD5 全通过 | 文件不可用 |
| 100K 清单覆盖 | 每个器官聚合后覆盖足够训练，目标 ≥80% | 数据太稀疏，需要改删失建模或缩小集合 |
| Liver 重建一致性 | 与公开 `Liver` 标签 Pearson r 建议 ≥0.80 | 分母、批次或处理流程有误 |
| 动物重复性 | 接近本审计表中的器官范围 | 重建方法未复现论文数据性质 |
| 官方测试性能 | Brain 约 0.61、spinal 约 0.65 为参考，不要求逐点相同 | 若明显更低，先修数据再换模型 |
| 数据泄漏 | 7-mer 不跨 train/test 重复；动物 4 不参与调参 | 否则结果不可答辩 |

若 Liver 一致性无法达到合理水平，主路线应暂停；届时只能做“包装 + 肝脏去靶向”的处理表项目，不能声称已经建立 CNS 序列模型。

## 11. 四周项目应如何调整

### 第 1 周：数据工程优先

- 固化 69-run manifest；
- 下载校验并流式计数；
- 先完成 Liver + Virus DNA，利用公开 Liver 列验证整条重建链；
- 通过后再下载脑、脊髓、心、肾。

### 第 2 周：单头模型与严格验证

- 包装：Ridge / Random Forest / 官方 LSTM 或等价复现；
- 器官：分别训练 5 个头；
- Animal 4 盲测，外加序列距离切分；
- 报 Pearson r、R²、MAE、bootstrap 区间和检测覆盖。

### 第 3 周：虚拟筛选

- 约束生成 7-mer；
- 包装硬门槛；
- 脑/脊髓奖励、肝/心/肾惩罚；
- 权重 ±0.10 敏感性；
- 帕累托、多样性和训练距离。

### 第 4 周：候选与答辩

- 24–30 条候选，分偏 CNS、偏低肝、折中三组；
- 输出逐候选预测、不确定性、最近训练序列和风险注释；
- 明确“早期小鼠器官 DNA 分布代理”的声明边界；
- 不声称替代 Zolgensma 或恢复 SMA 运动功能。

## 12. 本轮已经交付的可复现能力

仓库新增以下命令：

```bash
# 审计官方 GitHub release
aav9-sma audit-fit4function data/raw/fit4function_official \
  --output docs/audit_data/fit4function_release_audit.json

# 从 NCBI 实时生成 SRA manifest
aav9-sma fetch-sra-manifest \
  --bioproject PRJNA1131359 \
  --output-csv docs/audit_data/fit4function_sra_manifest.csv \
  --output-summary docs/audit_data/fit4function_sra_summary.json

# 对单个 FASTQ 做保守 21-nt / 7-mer pilot
aav9-sma count-fastq path/to/sample.fastq.gz \
  --whitelist-csv data/raw/fit4function_official/data/fit4function_library_screens.csv \
  --output docs/audit_data/sample_pilot.json

# 对公开 100K 带序列表做基线
aav9-sma benchmark-fit4function \
  data/raw/fit4function_official/data/fit4function_library_screens.csv \
  --models ridge random_forest \
  --output docs/audit_data/baseline_metrics.csv
```

当前 FASTQ 提取器是为了验证链路的轻量实现。完整重建阶段将增加与论文一致的 Bowtie2 路径、每 run 稀疏计数表、RPM/重复聚合、分母选择验证和动物留出数据集生成。

## 13. 一句话给队友或评委

> Fit4Function 的公开处理表足以直接训练包装和肝相关模型，却隐藏了多器官标签与 7-mer 的映射；我们通过审计 270 个公开 SRA runs，并实跑脊髓 FASTQ，确认可从原始 reads 重建这层映射。因此本项目的第一项真正工作不是“换一个更复杂的神经网络”，而是建立一条可验证的多器官序列标签重建管线，再在严格动物留出下做 CNS 奖励、肝脏和其他器官惩罚的多目标设计。

## Sources

[^1]: Eid, F.-E. et al. “Systematic multi-trait AAV capsid engineering for efficient gene delivery.” *Nature Communications* 15, 6602 (2024). [PMC full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11297966/) · [DOI](https://doi.org/10.1038/s41467-024-50555-y)
[^2]: Deverman Lab / vector-engineering. [Official Fit4Function GitHub repository](https://github.com/vector-engineering/fit4function).
[^3]: Chen, A. T. [vector-engineering/fit4function: v1.0.0](https://zenodo.org/records/8401253). Zenodo. DOI: 10.5281/zenodo.8401253.
[^4]: Eid, F.-E. et al. [Fit4Function source data files](https://zenodo.org/records/8388031). Zenodo. DOI: 10.5281/zenodo.8388031.
[^5]: NCBI. [BioProject PRJNA1131359](https://www.ncbi.nlm.nih.gov/bioproject/1131359), “Systematic Multi-Trait AAV Capsid Engineering for Efficient Gene Delivery.”
[^6]: Eid, F.-E. et al. [Supplementary Information](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-024-50555-y/MediaObjects/41467_2024_50555_MOESM1_ESM.pdf).
[^7]: NCBI SRA. [SRR29692586](https://www.ncbi.nlm.nih.gov/sra/SRR29692586), Hammerhead spinal-cord biodistribution technical replicate.
[^8]: vector-engineering/fit4function. [BSD 3-Clause license](https://github.com/vector-engineering/fit4function/blob/main/LICENSE).

# 运动神经元靶向模块：数据、模型与验证方案

> 本文整理项目中关于“从脊髓组织富集进一步走向脊髓运动神经元转导”的讨论，作为后续 dry-lab 实施和实验合作的技术交接文档。

## 执行摘要

当前最合理的补充版本不是立即声称已经能够预测运动神经元转导，而是在现有器官级模型之外建立一个独立、可升级的证据模块：

1. 保留现有 `spinal_cord_mouse` 作为组织级富集代理；
2. 使用公开单细胞/多组学数据建立运动神经元参考分类器；
3. 为每条外部证据记录 capsid、物种、路线、剂量、载荷、时间点和证据等级；
4. 只有拿到“capsid 序列 + 运动神经元标签 + 细胞级转导读数”后，才训练真正的 motor-neuron head；
5. 使用 study-level 外部留出验证，避免把同一动物或同一研究中的细胞随机拆分造成数据泄漏。

建议将这一补充版本命名为 **Motor-neuron Evidence Extension v0.1**。v0.1 的交付目标是数据契约、文献/数据清单、参考分类器设计和模型接口，不把间接证据包装成已经验证的运动神经元嗜性。

## 1. 研究问题与证据边界

当前项目的问题是：

> 在保留 AAV9 (K449R) 7-mer 插入衣壳包装潜力的前提下，能否提高脑/脊髓富集代理，同时降低肝脏和其他器官的富集代理？

这不是一个已经完成的“运动神经元靶向”模型。当前 Fit4Function 标签主要是小鼠器官层面的 vector-genome biodistribution proxy，测量的是整块器官中 AAV 相关 DNA 的相对丰度，而不是某种细胞中的功能表达。

必须区分：

~~~text
组织到达   = AAV 在脑或脊髓中出现
细胞进入   = AAV 进入运动神经元
细胞表达   = 载荷在运动神经元中表达
细胞特异性 = 相对于其他细胞，运动神经元中表达更高
功能效果   = 运动神经元功能或疾病表型被改变
~~~

项目当前主要处于第一层。运动神经元模块的目标是逐步增加第二至第四层证据，而不是把 spinal-cord 标签直接改名为 motor-neuron 标签。

### 当前模型可以说什么

- 某个 7-mer 在现有数据分布下具有较高的脑/脊髓富集预测。
- 某个候选的肝脏富集预测较低。
- 候选在当前序列距离、重复实验和 Animal 4 留出策略下具有一定的计算稳定性。
- 脑和脊髓是两个独立预测头，不能因为二者都属于 CNS 就默认它们等价。

### 当前模型不能说什么

- “预测进入运动神经元”。
- “运动神经元特异性衣壳”。
- “可用于 SMA 的运动神经元靶向治疗”。
- “在人或猕猴中具有与小鼠相同的运动神经元嗜性”。
- “低肝脏预测等于低肝毒性”。

当前脑/脊髓 Pearson 约为中等水平，且动物间重复性本身只有中等上限。因此模型适合做候选优先级排序，不足以支持细胞类型级结论。

## 2. 什么才算真实的运动神经元数据

| 级别 | 数据类型 | 能回答的问题 | 当前项目定位 |
| --- | --- | --- | --- |
| L1 | 整块脊髓 vector genome DNA | AAV 是否到达脊髓 | 现有 Fit4Function 主体 |
| L2 | 整块脊髓 RNA、蛋白或 reporter | 脊髓中是否有表达 | 外部参考或后续实验 |
| L3 | 细胞类型共定位、单细胞/单核表达、vector barcode | 是否进入运动神经元，哪些细胞表达 | 运动神经元模型的最低直接证据 |
| L4 | 运动神经元功能、疾病表型、独立组织学验证 | 进入是否产生功能效果 | SMA/ALS 验证阶段 |

“真实运动神经元数据”至少应包含运动神经元身份和该细胞中的转导或表达读数。只有整块脊髓 qPCR 或 bulk DNA 时，不能重命名为 motor-neuron data。

## 3. 公开论文和数据入口

### 3.1 最直接的参考：Kussick et al., 2025

**Enhancer AAVs for targeting spinal motor neurons and descending motor pathways in rodents and macaque**

[PubMed](https://pubmed.ncbi.nlm.nih.gov/40403722/) · DOI: 10.1016/j.celrep.2025.115730

该工作：

- 生成了小鼠和猕猴脊髓的 single-cell multiome 数据；
- 定义了脊髓运动神经元及其亚型；
- 从细胞类型数据中发现候选 enhancer；
- 构建 enhancer AAV reporter 工具；
- 在啮齿类和猕猴中进行了运动神经元相关验证。

这是目前最接近“细胞类型分辨的脊髓运动神经元数据”的公开工作。它适合用于：

1. 建立运动神经元和亚型的参考标签；
2. 训练运动神经元分类器；
3. 定义 promoter/enhancer 与细胞类型的关系；
4. 为后续模型提供外部验证集。

重要限制：该工作核心是 enhancer AAV，不一定提供当前项目 AAV9 7-mer 变体与逐细胞转导读数的配对关系。因此除非获得完整的 capsid/construct 序列和 cell-level readout，否则不能直接训练当前的 7-mer -> motor-neuron transduction 模型。

### 3.2 AAV9 与脊髓运动系统

**Bravo-Hernandez et al., 2020, Nature Medicine**

[Spinal subpial delivery of AAV9 enables widespread gene silencing and blocks motoneuron degeneration in ALS](https://pubmed.ncbi.nlm.nih.gov/31873312/)

DOI: 10.1038/s41591-019-0674-1

可用于参考 AAV9 在猪和非人灵长类脊髓运动系统中的表达、组织分布和 ALS 功能终点。限制是给药方式为 subpial delivery，不能直接与项目的系统给药/器官富集模型合并训练。

**Rashnonejad et al., 2019, Molecular Therapy**

[Fetal Gene Therapy Using a Single Injection of Recombinant AAV9 Rescued SMA Phenotype in Mice](https://pubmed.ncbi.nlm.nih.gov/31543414/)

DOI: 10.1016/j.ymthe.2019.08.017

可用于定义 SMA 中的 SMN 表达、运动功能和疾病表型终点，但不是 sequence-linked 的 capsid library 数据。

### 3.3 其他路线和 capsid 参考

- **Hordeaux et al., 2015, Gene Therapy**：AAVrh10 的静脉/鞘内 CNS 递送和脊髓/DRG 分布，适合做路线对照。 [PubMed](https://pubmed.ncbi.nlm.nih.gov/25588740/)
- **Deverman et al., 2016, Nature Biotechnology**：AAV-PHP.B 系统性 CNS 递送基准，适合做 CNS capsid benchmark，不是运动神经元特异性数据。
- **Chan et al., 2017, Nature Neuroscience**：AAV-PHP.eB 的非侵入性 CNS 递送。必须考虑小鼠品系和 LY6A 依赖性，不能直接外推到人或猕猴。

### 3.4 数据库和数据获取路径

按以下顺序检查：

1. 论文正文的 Data Availability 和 supplementary tables；
2. GEO、SRA、ENA 的原始单细胞/单核数据；
3. NeMO Archive；
4. Allen Brain Cell Types、BICAN、CELLxGENE 的脊髓细胞类型图谱；
5. 论文作者的补充数据或合作共享。

建议检索：

~~~text
spinal motor neuron multiome
spinal motor neuron enhancer AAV
AAV barcode spinal cord single cell
AAV motor neuron single nucleus RNA
capsid variant cell type transduction spinal cord
~~~

不能因为论文写了 single-cell 就假设它包含可训练的 capsid 数据。下载前必须确认是否有 variant/capsid ID、序列和细胞级转导读数。

## 4. capsid -> motor neuron 模块如何接入当前系统

### 4.1 推荐的多任务结构

保留当前 capsid encoder 和组织级预测头，新增细胞类型级 head：

~~~text
7-mer / full capsid sequence
             |
             v
       shared capsid encoder
             |
     +-------+--------+------------+
     v       v        v            v
 packaging brain   spinal cord  liver/off-target
                                 |
                                 v
                     motor-neuron entry
                     motor-neuron expression
                     glia/off-target
~~~

运动神经元 head 还应接收实验背景变量：

~~~text
species
mouse strain
route
dose
timepoint
payload
promoter/enhancer
~~~

概念上目标是：

~~~text
P(motor-neuron transduction)
= f(capsid sequence, species, route, dose, timepoint, payload)
~~~

### 4.2 为什么不能只用现有脊髓标签

脊髓 bulk readout 同时受到运动神经元、感觉神经元、中间神经元、胶质细胞、血管细胞、血管内残留、细胞组成、路线和时间点影响。因此：

~~~text
spinal cord enrichment != motor-neuron transduction
~~~

不能将 spinal_cord_mouse 自动转换成 motor_neuron 标签，也不能通过 bulk Pearson 反推细胞类型级 Pearson。

## 5. 数据契约：什么数据才能训练这个 head

建议新增一张 cell-type evidence table，字段至少包括：

~~~text
study_id
animal_id
variant_id
capsid_sequence
peptide_7mer
insertion_site
species
strain
route
dose
timepoint
payload
promoter_or_enhancer
tissue
cell_type
cell_subtype
readout_type
readout_value
replicate
batch
source_accession
evidence_level
~~~

真正可以训练 capsid -> motor neuron 的最小组合是：

~~~text
capsid/variant sequence
+
motor-neuron cell label
+
该细胞中的 vector barcode、RNA、蛋白或 reporter readout
~~~

如果缺少 capsid 序列但有细胞类型数据，只能训练 motor-neuron classifier。

如果有 capsid 序列但只有整块脊髓数据，只能训练 spinal-cord proxy。

如果两者来自不同实验且没有统一 route、dose、payload 和时间点，不应直接连接成监督标签。

## 6. 没有直接数据时能做什么

### 6.1 运动神经元参考分类器

先使用公开脊髓单细胞图谱训练 motor neuron、sensory neuron、interneuron、astrocyte、oligodendrocyte、microglia 和 endothelial/perivascular cell 分类器。

常见的运动神经元相关标记包括 CHAT、SLC18A3/VAChT、MNX1/HB9、ISL1/2 等，但必须以物种、数据集和细胞亚型注释为准。

### 6.2 组织级到细胞级的弱监督

对于只有 tissue readout、但论文报告了运动神经元共定位比例或 reporter 图像的研究，可使用弱标签：

~~~text
high-confidence motor-neuron evidence
moderate motor-neuron evidence
tissue-only evidence
unknown
~~~

输出应命名为 motor_neuron_evidence_score，不能写成已经验证的 tropism。

### 6.3 迁移学习和多任务学习

使用 Fit4Function 训练 shared capsid encoder，再用少量真实细胞类型数据微调 motor-neuron head。组织级任务作为辅助任务，不作为运动神经元标签。

### 6.4 生物学先验

Kussick 等 enhancer AAV 数据可以作为细胞类型和调控元件先验；AAV-PHP.B/PHP.eB 可以作为 CNS 递送先验；AAV9/SMA/ALS 论文可以作为功能终点先验。这些先验必须和直接监督标签分开存储。

### 6.5 不应做的推断

~~~text
脊髓富集高
-> 假设运动神经元富集高
-> 生成 motor-neuron label
~~~

这会把模型从“有局限的组织级预测”变成“没有实验依据的细胞级断言”。

## 7. “官方数据”是否必须

不需要政府或官方机构的数据，但必须满足“可追溯、可复核、实验定义清楚”：

| 数据来源 | 是否可用 | 说明 |
| --- | --- | --- |
| 已发表论文的原始矩阵/补充表 | 可以 | 记录 PMID/DOI、版本和下载地址 |
| GEO/SRA/ENA/NeMO/CELLxGENE | 可以 | 保留 accession 和处理脚本 |
| 论文作者共享的未发表数据 | 可以 | 明确使用许可和数据字典 |
| 只有论文图片、没有数值 | 只能作低等级证据 | 不应直接当精确监督标签 |
| 只有整块脊髓 qPCR | 只能作 tissue proxy | 不能称为 motor-neuron label |
| 模拟数据/规则生成标签 | 不能作真实生物学标签 | 只能测试代码和模型接口 |

如果目标是发表 capsid 对运动神经元转导的预测，最好还要有独立研究或实验数据作外部验证。仅用同一篇论文的数据训练和验证，容易产生研究级泄漏。

## 8. 验证方案

### 8.1 数据切分

优先使用 study-level leave-one-study-out：

~~~text
训练：若干研究、若干路线
验证：完全未参与训练的研究
测试：不同研究/物种/路线的外部数据
~~~

不能只随机切分细胞，因为同一动物、同一批次的细胞会同时出现在训练集和测试集，导致结果过于乐观。

### 8.2 必须报告的指标

- cell-level AUROC/AUPRC：识别运动神经元的能力；
- variant-level Pearson/Spearman：capsid 预测与运动神经元 readout 的相关性；
- calibration：预测概率是否可信；
- study-level performance：跨论文、跨批次稳定性；
- species/route shift：跨小鼠、猕猴和给药路线的下降幅度；
- negative-control specificity：胶质细胞、肝脏和其他 off-target 是否被误判。

### 8.3 建议的输出名称

在直接细胞级标签出现前，使用：

~~~text
pred_spinal_cord_proxy
motor_neuron_reference_score
motor_neuron_evidence_score
motor_neuron_head_uncertainty
domain_shift_warning
~~~

只有在存在 capsid 序列、细胞类型标签和细胞级转导读数，并且通过外部研究验证后，才可使用 pred_motor_neuron_transduction 和 pred_motor_neuron_specificity。

## 9. SMA 和其他疾病的应用边界

### SMA

运动神经元模块最自然的应用是 SMA。以下内容必须分开：

~~~text
capsid delivery
SMN1/SMN expression
motor-neuron survival
neuromuscular junction
motor behavior
~~~

不要将 capsid 富集预测直接等同于 SMA 疗效。

### ALS

ALS 的 SOD1、TDP-43、C9orf72 等模型可以提供运动神经元退行、存活和功能终点，但病理状态会改变细胞组成和血脑屏障。ALS 数据适合做疾病状态外部验证，不应直接与健康小鼠 bulk biodistribution 合并为同一标签。

### 其他可考虑疾病

- 脊髓损伤：运动通路重建和细胞类型特异表达；
- 遗传性痉挛性截瘫：长程运动通路和脊髓神经元；
- 部分运动神经元病和神经肌肉接头疾病。

疾病数据的作用是测试模块在病理条件下是否保持合理性，而不是单纯扩大训练集数量。

## 10. Dry-lab 实施路线

### Phase 0：数据审计

- [ ] 建立 cell_type_evidence 数据契约。
- [ ] 给每条外部记录记录 PMID/DOI、accession 和许可信息。
- [ ] 标记 capsid 序列是否真实可得。
- [ ] 标记 readout 是 DNA、RNA、蛋白、reporter 还是功能数据。
- [ ] 标记 route、dose、timepoint、payload 和 species。

### Phase 1：参考模型

- [ ] 下载并整理 Kussick 2025 相关数据。
- [ ] 建立运动神经元和非运动神经元参考分类器。
- [ ] 在 Allen/BICAN/CELLxGENE/NeMO 数据上做外部验证。
- [ ] 将结果保存为 reference evidence，而不是当前 capsid target。

### Phase 2：capsid head

- [ ] 复用现有 7-mer encoder。
- [ ] 保留 packaging、brain、spinal cord、liver 和 off-target heads。
- [ ] 只有在获得 sequence-linked cell-level 数据后新增 motor-neuron head。
- [ ] 加入 species、route、dose、timepoint、payload 的 context encoder。
- [ ] 使用 study-level split 和 uncertainty calibration。

### Phase 3：外部验证

- [ ] 选择完全未参与训练的研究作为测试集。
- [ ] 至少包含一个不同批次或不同实验室的数据源。
- [ ] 尽量包含一个不同物种或不同给药路线的数据源。
- [ ] 对 motor-neuron score、spinal-cord proxy 和 off-target score 分别报告结果。

### Phase 4：实验合作或验证

- [ ] 与有 AAV barcode、单细胞测序或组织学共定位能力的实验室合作。
- [ ] 预先固定细胞类型标记、主要 readout、阴性对照和 stop/go 规则。
- [ ] 先验证少量高、中、低预测候选，而不是只测 top hit。

## 11. 最终决策规则

在当前阶段，候选只能按以下层级表述：

~~~text
层级 1：spinal-cord enrichment candidate
层级 2：motor-neuron biological-prior candidate
层级 3：motor-neuron evidence-supported candidate
层级 4：experimentally validated motor-neuron targeting candidate
~~~

当前项目可以稳定支持层级 1，并可通过公开单细胞数据建立层级 2。层级 3 需要 sequence-linked 的细胞类型证据；层级 4 需要独立实验验证。

最重要的原则是：**可以没有“官方数据”，但不能没有可追溯的实验数据；可以先做 dry-lab 模块，但必须明确它是先验、代理还是经过细胞级验证的预测。**

## 12. 与现有仓库文档的关系

- [项目定义](PROJECT_SPEC.md)：当前研究问题、组织级标签和声明边界。
- [数据契约](DATA_CONTRACT.md)：Fit4Function 字段、标签含义和审计要求。
- [Fit4Function 审计](FIT4FUNCTION_AUDIT.md)：现有数据覆盖和复现结果。
- [实验验证 SOP](EXPERIMENTAL_VALIDATION_SOP.zh-CN.md)：候选验证和 go/hold/stop 规则。
- [三路线共识](THREE_ROUTE_CONSENSUS.md)：不同数据来源和分析路线的交叉检查。

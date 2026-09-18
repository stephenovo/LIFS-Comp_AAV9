# Fit4Function 数据运行手册

这份文件只讲“怎么跑”。科学解释和审计结论见 [`FIT4FUNCTION_AUDIT.md`](FIT4FUNCTION_AUDIT.md)。

## 1. 安装项目

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

## 2. 获取官方处理数据

```bash
git clone --depth 1 \
  https://github.com/vector-engineering/fit4function.git \
  data/raw/fit4function_official
git -C data/raw/fit4function_official checkout 6bfc2ebfe4abcd45cc6fe737e4700242a5090fee
```

The commit is pinned in [`source_manifest.json`](source_manifest.json). Do not
silently replace it with the upstream moving `main` branch when reproducing a
published result.

`data/raw/` 已被 Git 忽略，原始数据不会误推到 GitHub。

## 3. 审计官方 release

```bash
aav9-sma audit-fit4function \
  data/raw/fit4function_official \
  --output docs/audit_data/fit4function_release_audit.json
```

这个命令检查：

- 每个序列表的行数、列、哈希和 7-mer 合法性；
- 缺失值和无穷值；
- modeling / assessment / Fit4Function 序列交集；
- 多器官 workbook 是否真的含 AA 序列；
- 动物重复性和公开模型预测相关性；
- 包装双高斯分布的审计边界。

## 4. 生成最新 SRA manifest

```bash
aav9-sma fetch-sra-manifest \
  --bioproject PRJNA1131359 \
  --output-csv docs/audit_data/fit4function_sra_manifest.csv \
  --output-summary docs/audit_data/fit4function_sra_summary.json
```

manifest 是下载和批次映射的唯一入口。不要手工复制 69 个 accession。
仓库已同时保存筛选后的
[`audit_data/fit4function_sra_target_manifest.csv`](audit_data/fit4function_sra_target_manifest.csv)，
其中只有本项目所需的 60 个器官 runs 和 9 个病毒库分母 runs。

先把 ENA 的 FASTQ 地址、压缩大小和 MD5 加入 manifest：

```bash
aav9-sma resolve-ena-fastq \
  docs/audit_data/fit4function_sra_target_manifest.csv \
  --output docs/audit_data/fit4function_ena_fastq_manifest.csv
```

第一轮只下载 12 个 Liver 和 prod2/prod3 的 6 个 Virus DNA runs，共约 6.26 GB：

```bash
aav9-sma download-ena-fastq \
  docs/audit_data/fit4function_ena_fastq_manifest.csv \
  --output-dir data/raw/sra_samples \
  --workers 6 \
  --project-roles Liver "Virus reference" \
  --exclude-alias-regex prod1 \
  --output data/interim/fit4function_liver_downloads.json
```

下载器支持 HTTP 断点续传，并对每个完成文件同时核对字节数与 ENA MD5。

若 ENA 节点限速，可使用仓库中的 NCBI SRA 流水线。先按 NCBI 官方说明下载
macOS/Linux 对应的 SRA Toolkit，并把 `--sra-bin` 指向其 `bin/`。脚本会对每个归档运行
`vdb-validate`，随后转换、Bowtie2 计数，并只在成功写出 counts 和 QC 后删除临时文件：

```bash
PYTHONPATH=src python scripts/process_sra_to_counts.py \
  docs/audit_data/fit4function_ena_fastq_manifest.csv \
  --roles Brain "Spinal cord" Heart Kidney \
  --sra-bin /path/to/sratoolkit/bin \
  --s3-ip CURRENT_S3_IP \
  --work-dir tmp/sra_pipeline \
  --index-prefix data/interim/bowtie2/fit4function \
  --whitelist-csv data/raw/fit4function_official/data/fit4function_library_screens.csv \
  --counts-dir data/interim/bowtie2_counts \
  --summaries-dir data/interim/bowtie2_summaries \
  --workers 3 \
  --threads-per-worker 2
```

`--s3-ip` 是网络 DNS 失效时的恢复参数，不应长期写死；每次运行前应重新解析官方
`sra-pub-run-odp.s3.amazonaws.com`。直接 IP 路径仍以 `vdb-validate` 作为完整性硬门槛。
同一脊髓 run 的 ENA 与 SRA 路径已经做过逐 7-mer 交叉验证，100,000 行计数完全一致；
机器结果见 [`audit_data/SRR29692586_ena_sra_crosscheck.json`](audit_data/SRR29692586_ena_sra_crosscheck.json)。

## 5. 单个 FASTQ pilot

下面是本次已验证的脊髓 run。FASTQ 可以从 ENA 镜像按 accession 获取：

```bash
curl -L \
  https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR296/086/SRR29692586/SRR29692586.fastq.gz \
  -o data/raw/SRR29692586.fastq.gz

aav9-sma count-fastq \
  data/raw/SRR29692586.fastq.gz \
  --whitelist-csv data/raw/fit4function_official/data/fit4function_library_screens.csv \
  --counts-output data/interim/SRR29692586_counts.csv \
  --output docs/audit_data/SRR29692586_spinal_cord_pilot.json
```

下载后必须核对 ENA 提供的 MD5。本 run 的预期 MD5 是：

```text
d9bbf63bf1c40b8d7cfdb17f7aa1ea3a
```

## 6. 跑处理表基线

```bash
aav9-sma benchmark-fit4function \
  data/raw/fit4function_official/data/fit4function_library_screens.csv \
  --models ridge random_forest \
  --output docs/audit_data/fit4function_baseline_metrics.csv
```

这是随机留出 sanity check，不是最终评估。最终版本需要序列距离切分和动物留出。

## 7. 论文参数 Bowtie2 与 Liver 重建

先安装 Bowtie2，然后建立论文 149 nt 短参考的索引：

```bash
brew install bowtie2

aav9-sma prepare-bowtie2-reference \
  --output-prefix data/interim/bowtie2/fit4function \
  --output docs/audit_data/fit4function_bowtie2_reference.json
```

批量比对 18 个第一轮 FASTQ。SAM 通过内存流直接解析，不会写入磁盘：

```bash
aav9-sma count-bowtie2-batch \
  docs/audit_data/fit4function_ena_fastq_manifest.csv \
  --fastq-dir data/raw/sra_samples \
  --index-prefix data/interim/bowtie2/fit4function \
  --whitelist-csv data/raw/fit4function_official/data/fit4function_library_screens.csv \
  --counts-dir data/interim/bowtie2_counts \
  --summaries-dir data/interim/bowtie2_summaries \
  --workers 4 \
  --threads-per-worker 2 \
  --project-roles Liver "Virus reference" \
  --exclude-alias-regex prod1
```

最后同时计算“全部合格 7-mer reads 分母”和“公开 100K 白名单 reads 分母”，并遍历
Animal 1–4、Animal 1–3、单动物以及 prod2/prod3 病毒分母，避免手工挑一个最顺眼的组合：

```bash
aav9-sma reconstruct-liver \
  docs/audit_data/fit4function_ena_fastq_manifest.csv \
  --counts-dir data/interim/bowtie2_counts \
  --summaries-dir data/interim/bowtie2_summaries \
  --public-screens data/raw/fit4function_official/data/fit4function_library_screens.csv \
  --exclude-alias-regex prod1 \
  --output-reconstruction data/processed/fit4function_liver_reconstructed.csv.gz \
  --output-metrics docs/audit_data/fit4function_liver_validation_metrics.csv \
  --output-qc docs/audit_data/fit4function_liver_run_qc.csv
```

## 8. 完整 raw-data 管线顺序

本轮 10 步均已完成：

1. 从 manifest 只选 `Hammerhead`；
2. 先下载 12 个 Liver runs + 9 个 Virus DNA runs；
3. 按论文参数 Bowtie2 比对；
4. 提取 21 nt、Q20 过滤并翻译；
5. 只保留公开 100K AA 清单；
6. 输出每 run 稀疏 counts；
7. RPM、技术重复平均、log2 enrichment；
8. 与公开 `Liver` 标签做一致性验收；
9. 通过后再扩展到 Brain、Spinal cord、Heart、Kidney；
10. Animal 1–3 训练，Animal 4 作为 development held-out test；它已在模型比较中被查看，不能称为最终盲测。

最终建模表使用 60 个器官 runs 与 prod2 的 3 个病毒库 runs。prod3 只参与 Liver
分母敏感性分析；prod1 在 prod2 已通过公开 Liver 验收后无需下载。

## 9. 重建 5 个器官的逐动物标签

```bash
aav9-sma reconstruct-multiorgan \
  docs/audit_data/fit4function_ena_fastq_manifest.csv \
  --counts-dir data/interim/bowtie2_counts \
  --summaries-dir data/interim/bowtie2_summaries \
  --virus-round 2 \
  --output-reconstruction data/processed/fit4function_multiorgan_reconstructed.csv.gz \
  --output-metrics docs/audit_data/fit4function_multiorgan_replicate_metrics.csv \
  --output-qc docs/audit_data/fit4function_multiorgan_run_qc.csv
```

输出表有 100,000 行、123 列，包含两种 RPM 分母定义、四只动物、Animal 1–3 聚合、
四动物聚合及 5 个器官的 log2 enrichment。

## 10. 严格基线

```bash
aav9-sma benchmark-multiorgan \
  data/processed/fit4function_multiorgan_reconstructed.csv.gz \
  --models ridge random_forest \
  --output docs/audit_data/fit4function_multiorgan_baseline_metrics.csv
```

该命令用 Animal 1–3 标签训练、Animal 4 标签测试，并从训练集中排除与测试序列只有
一个氨基酸差异的 7-mer。不要用普通随机切分结果替代这张表。

## 11. 多任务集成与虚拟筛选

先评价筛选实际使用的五模型集成：

```bash
aav9-sma benchmark-multitask-ensemble \
  data/processed/fit4function_multiorgan_reconstructed.csv.gz \
  --ensemble-size 5 \
  --output docs/audit_data/fit4function_multitask_ensemble_metrics.csv
```

再运行 100 万条候选的完整漏斗。全量排序表较大，写入 Git 忽略的 `artifacts/`；
帕累托表、30 条清单和 JSON 摘要可以提交：

```bash
aav9-sma screen-virtual \
  data/raw/fit4function_official/data/fit4function_library_screens.csv \
  data/processed/fit4function_multiorgan_reconstructed.csv.gz \
  --pool-size 1000000 \
  --ensemble-size 5 \
  --max-iter 80 \
  --output-ranked artifacts/virtual_screen_ranked.csv.gz \
  --output-pareto docs/audit_data/virtual_screen_pareto.csv \
  --output-shortlist docs/audit_data/virtual_screen_shortlist.csv \
  --output-summary docs/audit_data/virtual_screen_summary.json
```

`pred_pack_lcb` 是在距离-2 校准集上估计的单侧 95% 下界。`uncertainty_*` 是五个
MLP 的模型分歧，不是校准置信区间。筛选逻辑和限制见
[`VIRTUAL_SCREEN_REPORT.md`](VIRTUAL_SCREEN_REPORT.md)。

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
```

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

## 7. 下一版完整 raw-data 管线

后续实现顺序固定为：

1. 从 manifest 只选 `Hammerhead`；
2. 先下载 12 个 Liver runs + 9 个 Virus DNA runs；
3. 按论文参数 Bowtie2 比对；
4. 提取 21 nt、Q20 过滤并翻译；
5. 只保留公开 100K AA 清单；
6. 输出每 run 稀疏 counts；
7. RPM、技术重复平均、log2 enrichment；
8. 与公开 `Liver` 标签做一致性验收；
9. 通过后再扩展到 Brain、Spinal cord、Heart、Kidney；
10. Animal 1–3 训练，Animal 4 盲测。

不要在第 8 步通过前批量训练 CNS 模型。数据处理链没有得到外部验证时，模型分数没有意义。

# AAV9–SMA Tropism Screening

面向脊髓性肌萎缩症（SMA）应用场景的 AAV9 衣壳组织特异性计算设计项目。

本项目计划在 AAV9 表面插入位点设计 7-mer 候选，以包装能力为硬门槛，利用公开数据训练多任务模型，奖励小鼠脑/脊髓富集代理，惩罚肝脏和其他器官脱靶，最终输出 24–30 条值得进一步实验验证的候选。

## 研究边界

- `F_pack`：可生产/包装预测，作为硬门槛；
- `F_CNS`：小鼠脑与脊髓标签构成的中枢富集代理；
- `F_liv`：小鼠肝脏标签，人肝细胞标签作为辅助注释；
- `F_off`：心、肾等其他器官脱靶代理；
- `R_imm`：抗体足迹或免疫相关风险，仅作为注释。

本项目不把小鼠器官信号称为“人运动神经元靶向”，不把低肝脏预测称为“已经降低肝毒性”，也不声称替代 Zolgensma。计算结果用于缩小后续实验搜索空间。

## 技术栈

- Python 3.11+
- NumPy、pandas：数据处理
- scikit-learn：Ridge、Random Forest 等基线模型
- LightGBM：可选增强基线
- PyTorch：可选共享编码器多任务模型
- pytest、ruff：测试和代码质量
- GitHub Actions：持续集成

复杂模型不是第一周的默认选择。项目先用可解释的基线验证数据是否含有可学习信号，再决定是否使用神经网络。

## 代码结构

```text
.
├── .github/workflows/ci.yml
├── configs/default.yaml
├── data/
│   ├── raw/          # 原始数据，不提交到 Git
│   ├── interim/      # 中间数据，不提交到 Git
│   └── processed/    # 建模表，不提交到 Git
├── docs/
│   ├── PROJECT_SPEC.md
│   ├── DATA_CONTRACT.md
│   └── TEAM_WORKFLOW.md
├── notebooks/        # 探索性分析；正式逻辑进入 src/
├── src/aav9_sma/
│   ├── data/audit.py
│   ├── features/encode.py
│   ├── models/baseline.py
│   ├── screening/pareto.py
│   ├── screening/score.py
│   └── cli.py
├── tests/
└── pyproject.toml
```

## 本地启动

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
pytest
```

可选模型依赖：

```bash
python -m pip install -e ".[dev,models]"
```

## 命令行骨架

审计已经映射到标准字段的数据：

```bash
aav9-sma audit-data data/processed/fit4function.csv \
  --output artifacts/data_audit.json
```

对模型预测结果执行包装过滤、综合评分和帕累托筛选：

```bash
aav9-sma rank-candidates artifacts/predictions.csv \
  --packaging-threshold 0.50 \
  --output artifacts/ranked_candidates.csv
```

`--packaging-threshold` 的正式数值必须由训练数据和验证结果确定；示例中的 `0.50` 不是预设生物学阈值。

## 默认展示分

在各任务标签标准化到可比较尺度后：

```text
S = 0.45 × F_CNS − 0.35 × F_liv − 0.20 × F_off
SI = F_CNS / (F_liv + ε)
```

综合分只用于展示和候选漏斗。科学交付同时保留偏中枢、偏低肝和折中三类帕累托候选，并进行权重敏感性分析。

## 当前里程碑

1. Fit4Function 数据可行性审计；
2. 建立标准数据映射和无泄漏切分；
3. 训练四个任务的简单基线；
4. 比较独立模型与共享编码器多任务模型；
5. 生成候选并完成包装硬过滤；
6. 综合评分、帕累托、多样性和不确定性筛选；
7. 输出候选卡、报告和答辩材料。


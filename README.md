# circaTE 0.1.0 — R 开发版

用配对 RNA-seq 与 Ribo-seq 计数检验**预先指定周期下的相对 TE 谐波变化**。第一版以 R 包实现，直接使用经过回归验证的 edgeR、limma 和可选 DESeq2 引擎；提供命令行入口，可在工作站或集群中运行。无需当前项目的历史目录，也不需要联网读取数据。

这是可安装的本地研究开发版，不是已经获得普遍 FDR 保证的正式统计方法发布。作者与维护者为 Phil Zhang（philzhang@tamu.edu）。开源许可证尚待确认；当前以私有 GitHub 仓库分发开发版。

## 安装

已验证环境为 R 4.4.2、edgeR 4.4.2、limma 3.62.2、DESeq2 1.46.0。其他版本可能改变数值，输出记录实际版本；改变版本后应重新做回归和科学校准。

```r
# 若依赖尚未安装；软件包安装不需要下载研究原始数据。
install.packages(c("BiocManager", "statmod", "jsonlite", "digest"))
BiocManager::install(c("edgeR", "limma"))
# 可选计数对照：
BiocManager::install("DESeq2")
```

```sh
R CMD INSTALL circaTE_0.1.0.tar.gz
```

GitHub 仓库：[zhangfei1947/circaTE](https://github.com/zhangfei1947/circaTE)（私有，需访问权限）。v0.1.0 作为开发预发布提供源码包；尚未发布到 CRAN/Bioconductor。

## 最小可运行示例

```r
library(circaTE)
d <- example_te_data()
f <- fit_te(d, method = "edger", period = 24)
f
head(f$results)
te_qc(d)$genes
plot_te_gene(f, "gene001")
write_te_results(f, "example_results") # 目标须不存在或为空
```

示例是 120 基因的人工接口演示，不能用它的检出数作为性能验证。

## 输入约定

- `rna`、`ribo` 为相同 gene × sample 顺序的原始整数计数矩阵。两者应基于一致且生物学适当的 feature 定义，例如 RNA CDS 与 Ribo CDS/P-site counts。
- 列名必须是相同的生物配对 ID；如果原始 RNA 和 Ribo 文库名不同，先依据可靠配对清单显式映射。函数不根据列位置猜配对，不自动交集删基因。
- `metadata` 至少包含 `sample_id`、数值 `time_hours`、`batch`，顺序与矩阵一致。`batch` 是沿时间连续跟踪的培养批次；同 batch/time 不能重复。每次拟合只处理一个生物组。
- 默认两种 assay 分别做 TMM，不自动进行丰度过滤，也不删除低计数观测。全部尝试的基因参与 BH。至少一个 assay 总计数为零的基因保留 p=1，并标为无信息。
- 若已筛选了矩阵，应提供从适当的全基因集合计算的自然对数 effective library size offsets；不要在小候选列表上重新估计全库归一化。输入原始 library size 当作 log offset 会导致错误，必须明确单位。

```r
d <- te_data(rna, ribo, metadata)
# 已有归一化：
d <- te_data(rna, ribo, metadata,
             offsets = list(rna = log_rna_library_sizes,
                            ribo = log_ribo_library_sizes))
# 读取 TSV，计数第一列为 gene_id，样本元信息字段如上：
d <- read_te_data("RNA.tsv", "Ribo.tsv", "samples.tsv")
```

## 方法与输出

| method | 用途与限制 |
|---|---|
| `edger` | 默认配对联合 edgeR v4 QL；吸收共同 RNA 变化，原始检验不显式建模序列相关 |
| `weighted_ar` | 观测精度加权＋培养内 AR1/nugget＋方差收缩；低计数 log 偏差和估计协方差不确定性仍在 |
| `unweighted_ar` | 匹配 H1 协方差估计的未加权 AR 对照，非历史 v2 H2 协方差版本 |
| `weighted_independent` | 保留精度加权但不拟合时间相关的对照 |
| `weighted_ar_ordinary` | 无方差收缩的加权 AR 对照 |
| `deseq2` | 原始配对 LRT，历史压力测试存在明显尾部误报，仅作研究对照 |
| `deseq2_finite` | 启发式有限自由度 F 参考，较保守，不是精确小样本检验 |

默认包含 TE batch 和线性趋势。`period=24` 是主分析；22 小时另作敏感性，不自动取两个周期或多个模型中最小 p。H2/adaptive 未进入本版正式 API：其非正弦检出收益伴随未解决的误报。

GLS AR 要求至少两个独立培养批次。加权模型的 categorical-time voom 精度估计还需要重复时间类别提供残差自由度；两批的实际时间若完全错开成不同类别，可能无法估权。软件仅检查实际使用的 H1 设计，修复了历史引擎多余的 H2 秩限制，因此支持 H1 可辨识的每 6 小时采样。此软件修复不修改历史冻结模型，原 4 小时设计有数值回归对照。公共入口还使用相对奇异值阈值 1e-8 拒绝数值混叠设计，例如 24 小时周期却只在 12 小时间隔采样，不能让浮点残余正弦列伪装成可辨识设计。非等距时间只允许非负 AR 相关；负相关限于规则时间格点。

`fit$results` 包含 `gene_id,p,q,cos,sin,amplitude_log2,phase_hours,informative,fit_failed,status`。振幅是 log2 TE 谐波半峰谷差；phase 是谐波峰相对时间零点的模周期时间，不包括趋势。q 在单次拟合的全部尝试基因内计算；跨组、跨方法和跨周期没有自动合并校正。

`fit$diagnostics` 保留收敛、相关参数边界、自由度等信息。未收敛不会偷偷切换模型；基因级失败保留 p=1。`te_qc()` 提供 RNA/Ribo 最小计数、低计数数量及培养内最长连续低谷；超过该培养最小采样间隔的缺口会中断连续长度；这是操作性观测计数定义，非低谷持续时长，不规则采样需结合实际时间解释。`plot_te_gene()` 画描述性观测轨迹，不把模型显著性当成生物学证明；其 `pseudocount` 仅改变图中的比值，不会重拟合，GLS 拟合仍使用固定 0.5 伪计数。

`write_te_results()` 输出逐基因结果和 QC TSV、完整 fit RDS、JSON 配置/版本/诊断、sessionInfo；拒绝覆盖非空目录。

## 命令行与 HPC

```sh
Rscript -e 'cat(system.file("scripts", "circaTE.R", package="circaTE"))'
# 用上面显示的实际路径替换下面占位路径：
Rscript /path/to/circaTE.R --rna RNA.tsv --ribo Ribo.tsv \
  --metadata samples.tsv --method edger --period 24 --out results_24h
```

CLI 默认按真实观测数据处理，不提供绕过模拟校准限制的开关。集群上应通过计算节点/作业调度器运行完整数据；本轮未修改 GRACE 环境或提交新远程作业。

## 实验性模拟校准

`build_te_calibrator()` 从独立全零模拟拟合构建映射；`calibrate_te()` 仅接受声明为 simulation 的数据，要求一致的模型、周期、趋势、时间设计、基因数、归一化、软件版本与 population context，并通过数据内容指纹和 dataset_id 双重检查拒绝重用训练数据。它无法验证用户提供的标签或全零声明是否真实，因此仍须正确设计模拟与独立留出。

```r
# null_fits 为研究者准备的独立全零模拟 fit 列表。
cal <- build_te_calibrator(null_fits, scenario = scenario_labels, all_null = TRUE)
heldout_calibrated <- calibrate_te(heldout_fit, cal, mode = "max")
```

不随软件捆绑 v2/v3 经验映射，不默认校准真实数据。不外推尾部，不输出 p=0；有限 p 分辨率可能限制 BH 单基因发现。即使模拟总体表现较好，也不保证低表达各子群的 FDR。

## 建模证据与软件验证

[三轮建模总结](inst/doc/MODELING_HISTORY.md) 给出 v1/v2/v3 的完整结果和推荐变化。**不同轮的模拟难度、过滤范围和波形不同，不能横向用 power 大小比较进步。** v3 中跨场景校准 edgeR 在预设压力分布下较稳健，但该结论不能包装成真实数据的通用统计保证。NCU02302 仍是对处理方式敏感的探索候选。

[适用范围](inst/doc/VALIDATION_SCOPE.md) 区分软件正确性、条件推断与科学验证。测试包括输入对齐、设计、零信息策略、历史数值回归、GLS 置换/不规则时间以及校准泄漏/并列概率。`validation/` 保存本地构建、安装、R CMD check 和集成验证记录，不随轻量源码包分发。

软件不包含原始研究数据，也不自动读取远程目录。原始三轮证据保留在本项目相邻目录，第一版不修改其冻结数值实现或结果。

## 获取 GitHub 源码

```sh
gh repo clone zhangfei1947/circaTE
R CMD INSTALL circaTE
```

仓库不包含研究原始计数、私人路径配置或大型本地验证产物。随仓库保留可运行测试、方法总结与 [验证摘要](VALIDATION.md)；部分与历史研究目录的数值对照测试在独立 checkout 中会明确跳过，核心合成数据测试仍可独立运行。

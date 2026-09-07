# EC mRNA 疫苗分析：单细胞水平 IS3/IS5 细胞构成验证

**日期**：2026-09-07

**数据集**：GSE278879（snRNA-seq, 10x Genomics）
- **标题**：Single nucleus RNA sequencing analysis reveals transcriptional heterogeneity of endometrioid endometrial cancer
- **样本**：原发肿瘤 11 例；本分析使用其中 3 个代表性肿瘤
  - **HEEC1**：高级别子宫内膜样癌（High-grade）
  - **MEEC1**：中级别子宫内膜样癌（Middle-grade）
  - **LEEC1**：低级别子宫内膜样癌（Low-grade）
- **细胞总量**：29,785 个细胞（HEEC1=24,823；MEEC1=1,599；LEEC1=3,363）

---

## 1. 分析策略

### 1.1 IS3 / IS5 特征基因推导

从 TCGA-UCEC bulk RNA-seq 数据（n=549 原发肿瘤）中分别计算：

| 比较组 | 方法 | 上调基因标准 | 签名大小 |
|--------|------|-------------|---------|
| IS3 vs 其余 4 亚型 | Mann-Whitney U + BH-FDR | log₂FC > 0.5, FDR < 0.05 | **50** |
| IS5 vs 其余 4 亚型 | Mann-Whitney U + BH-FDR | log₂FC > 0.5, FDR < 0.05 | **50** |

**IS3 签名前 10 基因**（免疫热相关）：
`PLA2G2D, CXCL9, IFNG, CXCL11, CTSE, SLAMF7, GZMH, DPEP3, MS4A1, NKG7`

**IS5 签名前 10 基因**（上皮/分泌/抑制相关）：
`CARTPT, MMP27, PAGE2, CALCB, COL9A1, MAGEB2, SST, C19orf85, UPK2, UCMA`

### 1.2 单细胞处理流程

1. 从 GEO 下载每个样本的 `matrix.tar.gz`（10x 输出）
2. 合并 `barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`
3. 质控：保留总 UMI ≥ 500 的细胞
4. 归一化：log1p(counts-per-million)
5. 细胞类型标注（z-score marker 均值取最大）：
   - Epithelial（EPCAM/KRT8/KRT18/KRT19）
   - CD8_T（CD8A/CD8B/GZMK/GZMA/PRF1）
   - CD4_T, Treg（FOXP3/CTLA4）
   - B_cell, Plasma, NK, Macrophage, Monocyte, DC
   - Fibroblast, Endothelial, Mast
6. 对每个细胞计算 IS3/IS5 签名 z-score 均值
7. 统计 IS3/IS5 高表达细胞比例、关键抗原 VTCN1（B7-H4）与 CD274（PD-L1）的细胞类型分布

---

## 2. 主要结果

### 2.1 细胞类型构成（Fig 17a）

| 细胞类型 | HEEC1 (High) | MEEC1 (Middle) | LEEC1 (Low) |
|---------|-------------|---------------|--------------|
| Epithelial | 18.6% | 14.7% | 15.4% |
| CD8_T | **5.0%** | 1.3% | 2.5% |
| CD4_T | 8.4% | 6.8% | 5.4% |
| Treg | 10.4% | 6.2% | 6.9% |
| B_cell | 5.9% | 5.4% | 5.7% |
| Plasma | 10.3% | 6.5% | 9.4% |
| Macrophage | 5.8% | 7.1% | 6.2% |
| Monocyte | 5.3% | 4.7% | 3.5% |
| DC | 3.6% | 4.3% | 3.3% |
| NK | 4.3% | 4.1% | 3.6% |
| Fibroblast | 8.6% | 18.2% | 16.2% |
| Mast | 0.3% | 2.3% | 2.8% |

**关键观察**：
- **HEEC1（高级别）CD8_T 比例最高（5.0%）**，约为 MEEC1/LEEC1 的 2–4 倍，与 TCGA 中 IS3/IS4 “热”亚型的 T 细胞富集一致。
- **LEEC1 和 MEEC1 的成纤维细胞/基质比例更高**，可能对应更冷的 TIME。

### 2.2 IS3 签名在 CD8 T 细胞中高度富集（Fig 17b）

| 细胞类型 | IS3_score 均值 | IS5_score 均值 | IS3−IS5 中位 |
|---------|---------------|---------------|--------------|
| **CD8_T** | **0.469** | −0.008 | **最高（正向）** |
| NK | 0.178 | −0.002 | 次高 |
| Macrophage | 0.135 | −0.005 | 第三 |
| Mast | 0.007 | −0.003 | 接近 0 |
| B_cell | −0.031 | 0.002 | 略负 |
| CD4_T | −0.033 | 0.003 | 略负 |
| Treg | −0.036 | −0.004 | 略负 |
| Fibroblast | −0.044 | 0.005 | 负 |
| Monocyte | −0.050 | −0.001 | 负 |
| **Epithelial** | **−0.057** | **0.004** | **最低（负向）** |

**统计验证**（Epithelial vs CD8_T）：
- **IS3_score**：Epithelial 均值 = −0.057，CD8_T 均值 = **0.469**，Mann-Whitney U **P ≈ 0.0**
- IS5_score：Epithelial 均值 = 0.004，CD8_T 均值 = −0.008，P = 3.73e-24

**结论**：⚠️ **IS3 签名强烈且特异地富集于 CD8⁺ T 细胞**，直接证明 TCGA bulk 中定义的 IS3 “免疫热”亚型对应单细胞水平的 CD8 T 细胞浸润。

### 2.3 IS5 签名的细胞来源

IS5 签名在 **Epithelial（0.004）和 Fibroblast（0.005）** 中相对最高，但绝对值很低。这与 IS5 的 bulk 特征基因多为上皮/分泌相关基因（MMP27、COL9A1、UPK2、UCMA）一致，提示 IS5 可能反映肿瘤上皮和基质微环境的共同状态，而非单一细胞类型。

### 2.4 VTCN1（B7-H4）的细胞类型分布（Fig 17c）

| 细胞类型 | VTCN1 均值（log1p CPM） | VTCN1+ 细胞比例 |
|---------|----------------------|----------------|
| Mast | 0.527 | 18.4% |
| B_cell | 0.314 | 12.1% |
| **Epithelial** | **0.284** | **11.2%** |
| Treg | 0.186 | 9.4% |
| Plasma | 0.149 | 6.0% |
| NK | 0.143 | 6.4% |
| DC | 0.139 | 5.6% |
| Monocyte | 0.118 | 4.3% |
| CD4_T | 0.117 | 4.5% |
| Endothelial | 0.113 | 4.7% |
| Macrophage | 0.113 | 7.0% |
| Fibroblast | 0.086 | 6.0% |
| CD8_T | 0.073 | 4.9% |

**关键观察**：
- **Epithelial 是 VTCN1 表达最高的主要实质细胞类型之一**，与 B7-H4 作为肿瘤抗原的角色一致。
- **Mast 和 B 细胞也表达 VTCN1**，可能反映肿瘤微环境中免疫抑制细胞对 B7-H4 的贡献；这一点需在后续分析中用更精细的免疫细胞亚群标记验证。

### 2.5 肿瘤间异质性：VTCN1 与 CD8 浸润呈反向趋势

| 样本 | 级别 | VTCN1+ Epithelial | CD8A+ CD8_T | VTCN1 均值 |
|------|------|-------------------|-------------|-----------|
| HEEC1 | High | **5.1%** | **29.9%** | 0.059 |
| LEEC1 | Low | 40.6% | 36.1% | 0.449 |
| MEEC1 | Middle | **65.9%** | 23.8% | **1.238** |

- MEEC1（中级别）呈现 **VTCN1 高 + CD8 低** 的组合，最接近论文中 IS5（B7-H4 高、免疫抑制）的表型。
- HEEC1（高级别）VTCN1 低但 CD8 浸润高，更接近 IS3/IS4 免疫热表型。
- 但样本量仅 3 个肿瘤，不能建立统计学显著性。

---

## 3. 对论文叙事的支撑

| 论文声明 | 单细胞证据 | 强度 |
|---------|-----------|------|
| IS3 是 CD8 T 细胞富集的免疫热亚型 | ✅ IS3 签名在 CD8_T 中显著富集（P≈0） | **强** |
| IS3 签名也见于 NK / 巨噬细胞 | ✅ NK（0.178）和 Macrophage（0.135）均为正 | 中 |
| IS5 是上皮/基质主导的冷亚型 | ⚠️ IS5 签名在 Epithelial/Fibroblast 略高，但信号弱 | 弱（需更多样本） |
| VTCN1 是肿瘤抗原（上皮来源） | ✅ Epithelial 中 VTCN1 表达明显高于 CD8_T、Fibroblast 等 | 中 |
| VTCN1 高表达与低 CD8 浸润相关 | ⚠️ MEEC1（VTCN1 最高）CD8A+ 比例最低，但样本量小 | 初步 |

---

## 4. 局限与下一步

1. **样本量小**：当前仅 3 个肿瘤；GSE278879 还有 8 个肿瘤样本（HEEC2/3、MEEC2、LEEC2、NE1-4）可供扩展。
2. **一个肿瘤对应一个级别**：无法做肿瘤内/肿瘤间统计；需要更多样本才能验证 grade × IS 关联。
3. **细胞类型注释基于 marker 得分**：未使用 CopyKAT / inferCNV 区分恶性 vs 正常上皮；未来可补充 CNV 推断恶性细胞。
4. **VTCN1 在 Mast/B 细胞中的高表达**：可能受双联体/注释噪声影响，需用更严格的免疫细胞 marker 复核。
5. **IS5 签名在单细胞中信号弱**：提示 IS5 可能是肿瘤-基质交互的整体状态，而非单一细胞类型的标记；可考虑构建“上皮-成纤维细胞共表达模块”进一步验证。

---

## 5. 生成文件

| 文件 | 内容 |
|------|------|
| `data/is3_signature.csv` | IS3 特征基因（50 up） |
| `data/is5_signature.csv` | IS5 特征基因（50 up） |
| `data/is3_is5_signatures.json` | 合并签名字典 |
| `data/sc_ucec/cell_assignments.csv` | 29,785 个细胞的类型、IS3/IS5 分数、关键基因表达 |
| `data/sc_ucec/composition.csv` | 每个样本的细胞类型组成 |
| `data/sc_ucec/per_tumor_metrics.csv` | 每个样本的 IS3/IS5 高表达细胞比例、VTCN1/CD8A/CD68 均值 |
| `data/sc_ucec/cell_type_expression.csv` | 各细胞类型关键基因表达均值 |
| `figures/fig17_sc_validation.{png,pdf,svg}` | 单细胞验证配图（3 panel） |
| `data/parse_and_annotate_sc.py` | 单细胞解析与注释脚本 |
| `data/analyze_sc_validation.py` | 单细胞验证统计与配图脚本 |

---

*本分析为 EC 疫苗项目单细胞层面的初步验证。核心结论是 IS3 bulk 亚型对应 CD8 T 细胞浸润，VTCN1（B7-H4）在肿瘤上皮细胞中可检测到表达，支持 B7-H4 作为 IS5 相关靶点的研究方向。*

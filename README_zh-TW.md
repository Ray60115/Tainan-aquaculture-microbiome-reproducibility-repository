# 台南養殖池微生物相分析重現性資料夾

這份資料夾整合了本研究正式投稿稿件所需的完整分析脈絡，包括：

- PacBio HiFi-16S workflow v0.9 的完整程式快照；
- 本研究實際使用的 Nextflow 設定、參數與 204 筆輸入的 QC；
- `204 → 190 → 168` 的樣本對照；
- `29,461 → 25,219 → 23,896` 的 ASV 篩選紀錄；
- Shannon diversity、rarefaction、Bray–Curtis PCoA、限制式置換檢定、
  PERMANOVA／PERMDISP、環境向量、turnover 與變異分割；
- random forest、六種演算法敏感度分析，以及排除 T1 的 T2–T4 分析；
- PICRUSt2 的輸入準備、推估流程、功能分析與敏感度分析；
- 主要圖件重建程式與可供核對的 reference results。

## 已確認的樣本與 ASV 流程

- workflow 共有 204 筆輸入。
- 扣除 8 筆 `-G`、1 筆 Water blank、3 筆 probiotic reference、2 筆
  seawater reference 後，正式生物樣本 library 為 190 筆。
- 其中 22 個 pond-month-layer profile 各有兩筆來源 library；平均後可精確
  重建 168 個 profile，所有 ASV cell 的最大差異為 0。
- 原始表有 29,461 個 ASV；4,242 個未出現在 168 個正式 profile，另排除
  1,323 個 mitochondrial／chloroplast／explicit eukaryotic ASV，最後保留
  23,896 個 prokaryotic ASV。
- workflow 雖列出 204 筆輸入，但 feature table 實際有 199 個樣本欄位；
  缺少的是 3 筆 probiotic reference 與 2 筆 seawater reference，均不屬於
  正式 190-library cohort。190 筆正式 library 全部存在。
- 正式 168-profile ASV 表只含這 168 欄；8 筆 `-G` 與 1 筆 Water control
  已另存，避免把非正式樣本誤看成分析樣本。

## 執行方式

安裝 Python 套件：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

重跑主要分類與統計分析：

```bash
bash run_core_analysis.sh
```

若只要快速核對稿件中的樣本數、ASV 數與封存結果，可執行：

```bash
python code/00_validate_manuscript_invariants.py
```

利用已封存的 PICRUSt2 primary outputs 重跑功能統計；若尚未先執行 core
analysis，程式會改用已封存的 reference turnover table：

```bash
bash run_functional_analysis.sh
```

若要完整重跑保守的 blank-screened 功能敏感度分析（PICRUSt2 pathway inference、
下游功能統計，以及 primary 與 screened 結果比較），請先啟用已記錄的 PICRUSt2
環境，再執行：

```bash
PICRUST2_PROCESSES=4 bash run_blank_screened_functional_sensitivity.sh
```

## 公開資料與完成的敏感度分析

BioProject 為 `PRJNA1513996`，SRA study 為 `SRP727697`。由 NCBI 回傳的
185 筆與後補 5 筆 metadata 已合併，190 筆正式 library 均有唯一的 SAMN
與 SRR，完整對照表已放在 `data/manifests/`。2026-10-02 已由未登入狀態
確認 BioProject 顯示完整 190 筆 SRA experiments，公開資料可正常查詢。

blank-screened PICRUSt2 sensitivity analysis 所需的 combined EC prediction
與重跑輸出均已補入。排除 109 個 blank-enriched ASVs 後保留 23,787 個 ASVs，
產生 555 條 MetaCyc pathways。與 primary 分析比較時，pathway Bray–Curtis
距離的 Spearman rho 為 0.999877（稿件精度為 `0.9999`），22 項可比較檢定的
0.05 顯著性分類均未改變，因此這項稿件敘述現在可由本資料夾獨立重現。

更完整的程式對照請見 `MANUSCRIPT_TO_CODE_MAP.md`，重現性限制與版本說明請見
`REPRODUCIBILITY_NOTES.md`，逐項核對結果請見
`MANUSCRIPT_REPOSITORY_AUDIT_20261001.md`。

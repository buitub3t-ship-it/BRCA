# OS–M2EFM replication: TCGA training, Kao/Terunuma validation

Bộ code này triển khai trực tiếp kiến trúc M2EFM cho **overall survival** bằng R hiện đại, không phụ thuộc package M2EFM cũ có `biocLite()`:

1. Map Terunuma Human Gene 1.0 ST transcript-cluster ID sang gene symbol và lấy median khi nhiều probe-set cùng gene.
2. Giữ gene chung giữa TCGA, Kao và Terunuma.
3. Dùng 550 CpG đầu của `PROBELIST.csv`.
4. Chạy MatrixEQTL trên TCGA expression–methylation matched tumors.
5. Chọn top 110 trans-m2eGenes theo `|beta|`, một association mạnh nhất trên mỗi gene.
6. Full internal model: trans expression + selected m2eQTL methylation → Cox-Ridge molecular risk → Cox model với stage + age.
7. External model: trans genes + cis replacement genes → Cox-Ridge molecular risk → Cox model với stage + age; validate trên Kao và Terunuma.
8. Lặp 100 random split 70/30 và fit final model bằng toàn bộ TCGA.

## Cấu trúc thư mục

Đặt folder code cạnh `csv_output`:

```text
project_parent/
├─ csv_output/
│  ├─ TCGA_BRCA_EXP.csv
│  ├─ TCGA_BRCA_METH.csv
│  ├─ TCGA_BRCA_CLIN.csv
│  ├─ KAO_BRCA_EXP.csv
│  ├─ KAO_BRCA_CLIN.csv
│  ├─ TERUNUMA_BRCA_EXP.csv
│  ├─ TERUNUMA_BRCA_CLIN.csv
│  └─ PROBELIST.csv
└─ M2EFM_OS_TCGA_KAO_TERU/
```

Hoặc đặt CSV ngay trong folder code. Hoặc đặt biến môi trường `M2EFM_DATA_DIR` trỏ tới folder chứa CSV.

## Chạy trong RStudio 4.5.x

Mở project folder làm working directory, sau đó:

```r
source("00_install_packages.R")  # chỉ chạy lần đầu
source("run_all.R")
```

Chạy từng phase khi cần debug:

```r
source("01_preflight_and_prepare.R")
source("02_identify_m2eqtls.R")
source("03_train_internal_full.R")
source("04_train_external_expression.R")
source("05_save_session_info.R")
```

## Output chính

- `results/01_data_audit.csv`: kiểm tra kích thước, ID overlap, complete OS.
- `results/02_trans_m2eqtl_selected.csv`: top trans associations.
- `results/02_expression_only_signature.csv`: signature dùng cho Kao/Terunuma.
- `results/03_internal_full_summary.csv`: TCGA Meth+Exp+Clin qua 100 split.
- `results/04_expression_external_summary.csv`: TCGA/Kao/Terunuma Exp+Clin qua 100 split.
- `results/04_final_expression_external_metrics.csv`: final model train toàn bộ TCGA.
- `results/04_final_expression_external_model.rds`: model, scaling và signature để predict lại.

## Trạng thái bộ dữ liệu hiện tại

- `PROBELIST.csv` có 7.318 probe; 7.313 có trong TCGA methylation.
- Không probe nào trong **top 550** bị thiếu, nên đủ để bắt đầu m2eQTL discovery.
- TCGA expression: 10.000 genes × 1.059 tumors.
- TCGA methylation: 7.313 probes × 878 samples; 748 tumor samples overlap expression, 743 có complete OS trong clinical.
- Kao: 327 samples, expression và clinical khớp hoàn toàn.
- Terunuma: 61 samples, expression và clinical khớp hoàn toàn; script sẽ map 33.297 transcript clusters về gene symbols.

## Giới hạn tái lập

Các CSV hiện tại là repo-style data, nhưng không hoàn toàn trùng cohort công bố: file có 10.000 genes và khoảng 1.039 TCGA expression samples có complete OS, trong khi Supplement mô tả 10.990 genes và 1.028 samples. Vì vậy code tái lập **workflow/model architecture** trên dữ liệu hiện có; C-index và exact signature có thể khác paper.

`APPLY_COMBAT_TO_CURRENT_MATRICES` mặc định là `FALSE`. TCGA và Kao hiện có gần như cùng gene-wise/global distribution, cho thấy đã được harmonize. Không nên chạy ComBat lần nữa trên hai matrix đã xử lý. Để tái lập preprocessing paper tuyệt đối, cần quay lại matrices trước ComBat của cả ba cohort và batch-correct lại cùng lúc.


## Smoke test trước khi chạy 100 splits

Để kiểm tra package, annotation và format dữ liệu, tạm đổi trong `config.R`:

```r
N_MONTE_CARLO_SPLITS <- 3L
```

Chạy `source("run_all.R")`. Khi toàn bộ phase hoàn tất, đổi lại `100L` và chạy lại phase 03–04.

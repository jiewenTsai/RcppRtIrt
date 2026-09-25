# PISA 2018 empirical illustration (SMI RT-IRT)

Files:
- `download_pisa2018.sh` — fetches and unzips `SPSS_STU_COG.zip` and `SPSS_STU_QQQ.zip` from
  webfs.oecd.org into `data/pisa2018/` (network access to that host is required; the data are not
  committed).
- `prepare_pisa2018.R [CNT] [domain] [min items]` — one country, one domain (default USA, math
  "CM"; reading is adaptive in 2018). Picks the most common complete block of administered and
  timed items, codes full credit as 1, converts item times to log seconds, and merges gender
  (`ST004D01T`), language at home (`ST022Q01TA`, 2 = other language) and `ESCS`. Variable names
  follow the PISA 2018 codebook and should be checked against it.
- `analyze_pisa2018.R <prepared.rds> [n_iter]` — SMI fits over eta in {0, .1, .25, .5, .75, 1} with
  gender in the theta and tau means; targets: female gap (in the model), other-language gap and ESCS
  slope (not in the model); risk-rule eta per target; tau-score screen by language and ESCS;
  item-level RT residual gaps by language (differential response time); leakage coefficient.
- `make_fake_pisa.R` — synthetic files with the same layout, used to test the pipeline.

Run:
```
bash pisa/download_pisa2018.sh
Rscript pisa/prepare_pisa2018.R USA CM 8
Rscript pisa/analyze_pisa2018.R data/pisa2018/prepared_USA_CM.rds 4000
```

Pipeline test on the fake files (n = 806, p = 10; multilingual students slower overall and on
items 1-3): the screen flags language (LM = 127) and not ESCS (p = .44), items 1-3 show residual
gaps with t ≈ 7-9, and the risk rule picks eta = 0.1 for the language gap, where the full model
moves the gap from -0.18 (cut) to -0.13.

#!/usr/bin/env bash
# download_pisa2018.sh — fetch the PISA 2018 student cognitive and questionnaire files (SPSS).
# Needs network access to webfs.oecd.org. The file names follow the OECD PISA 2018 database page
# (https://www.oecd.org/pisa/data/2018database/); if a link has moved, download the two SPSS zips
# by hand into data/pisa2018/ and rerun this script to unzip them.
set -euo pipefail
DIR="${PISA_DIR:-data/pisa2018}"
mkdir -p "$DIR"
for f in SPSS_STU_COG.zip SPSS_STU_QQQ.zip; do
  if [ ! -f "$DIR/$f" ]; then
    curl -fL --retry 3 -o "$DIR/$f" "https://webfs.oecd.org/pisa2018/$f"
  fi
  unzip -o -q "$DIR/$f" -d "$DIR"
done
ls -la "$DIR"

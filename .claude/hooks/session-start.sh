#!/bin/bash
# SessionStart hook for Claude Code on the web:
# installs R, Rcpp/RcppArmadillo and the pg (Polya-Gamma) headers so that
# sourceCpp("RtIrtGibbs.cpp") and the helper scripts work in cloud sessions.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

# CRAN is not reachable under the default network policy, so R packages come
# from Ubuntu's prebuilt r-cran-* binaries, and pg is built from GitHub.
APT_PKGS=(
  r-base-dev
  r-cran-rcpp
  r-cran-rcpparmadillo
  r-cran-mass
  r-cran-coda
  r-cran-testthat
  r-cran-lintr
)

missing=()
for p in "${APT_PKGS[@]}"; do
  dpkg -s "$p" >/dev/null 2>&1 || missing+=("$p")
done
if [ "${#missing[@]}" -gt 0 ]; then
  export DEBIAN_FRONTEND=noninteractive
  # Some preconfigured PPAs are blocked by the proxy; ignore their failures.
  apt-get update -qq || true
  apt-get install -y -qq --no-install-recommends "${missing[@]}"
fi

# pg (tmsalab/pg) is not on CRAN; pin the commit for reproducibility.
PG_REF="348b5d499ea75da9565488d64d532420802e84e2"  # pg 0.2.4
if ! Rscript -e 'quit(status = !requireNamespace("pg", quietly = TRUE))' >/dev/null 2>&1; then
  tmp="$(mktemp -d)"
  git clone --quiet https://github.com/tmsalab/pg.git "$tmp/pg"
  git -C "$tmp/pg" checkout --quiet "$PG_REF"
  R CMD INSTALL --no-test-load "$tmp/pg"
  rm -rf "$tmp"
fi

Rscript -e 'for (p in c("Rcpp", "RcppArmadillo", "pg", "MASS", "coda")) stopifnot(requireNamespace(p, quietly = TRUE))'

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

# apt's Rcpp (1.0.12) caps List::create() at 20 elements, which
# RtIrtGibbs.cpp exceeds, so build a newer release from GitHub.
RCPP_VERSION="1.1.2"
RCPP_REF="93987a4014cd6670efa4d25dd279b0648dc8b25f"  # tag 1.1.2
if ! Rscript -e "quit(status = packageVersion('Rcpp') < '$RCPP_VERSION')" >/dev/null 2>&1; then
  tmp="$(mktemp -d)"
  git clone --quiet --filter=blob:none https://github.com/RcppCore/Rcpp.git "$tmp/Rcpp"
  git -C "$tmp/Rcpp" checkout --quiet "$RCPP_REF"
  R CMD INSTALL --no-docs "$tmp/Rcpp"
  rm -rf "$tmp"
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

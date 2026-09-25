# prepare_pisa2018.R — one country, one domain: a complete block of scored responses and item times.
#
# Inputs (OECD PISA 2018 database, SPSS format), placed in data/pisa2018/:
#   CY07_MSU_STU_COG.sav  cognitive item data: scored response <item>S and total time <item>T
#   CY07_MSU_STU_QQQ.sav  student questionnaire: ST004D01T (gender, 1 = female), ST022Q01TA
#                         (language at home, 1 = language of the test), ESCS
# The names follow the PISA 2018 codebook; check them against the codebook shipped with the data.
# Reading is adaptive in 2018 (multistage), so the default domain is mathematics ("CM"), whose
# clusters are fixed forms.
#
# Output: data/pisa2018/prepared_<CNT>_<domain>.rds with Y (n x p, 0/1 full credit), logT (log
# seconds), female, ml (1 = other language at home), escs, items, and a selection log.
#
# Usage: Rscript pisa/prepare_pisa2018.R [CNT] [domain prefix] [min items]
suppressMessages({ library(haven) })
args <- commandArgs(trailingOnly = TRUE)
cnt <- if (length(args) >= 1) args[1] else "USA"
dom <- if (length(args) >= 2) args[2] else "CM"
min_items <- if (length(args) >= 3) as.integer(args[3]) else 8
dir <- Sys.getenv("PISA_DIR", "data/pisa2018")
cog_file <- file.path(dir, "CY07_MSU_STU_COG.sav")
qqq_file <- file.path(dir, "CY07_MSU_STU_QQQ.sav")
stopifnot(file.exists(cog_file), file.exists(qqq_file))

pat <- paste0("^", dom, "[0-9].*[ST]$")
cog <- read_sav(cog_file, col_select = c("CNT", "CNTSTUID", tidyselect::matches(pat)))
cog <- cog[as.character(cog$CNT) == cnt, ]
cat("students in", cnt, ":", nrow(cog), "\n")
s_cols <- grep(paste0("^", dom, "[0-9].*S$"), names(cog), value = TRUE)
stems <- sub("S$", "", s_cols)
stems <- stems[paste0(stems, "T") %in% names(cog)]
cat(dom, "items with both a score and a time:", length(stems), "\n")
S <- sapply(stems, function(s) as.numeric(zap_labels(cog[[paste0(s, "S")]])))
Tm <- sapply(stems, function(s) as.numeric(zap_labels(cog[[paste0(s, "T")]])))
ok <- !is.na(S) & !is.na(Tm) & Tm > 0

# the most common set of administered-and-timed items (a cluster or cluster pair)
key <- apply(ok, 1, function(r) paste(which(r), collapse = ","))
tab <- sort(table(key[nchar(key) > 0]), decreasing = TRUE)
sizes <- sapply(strsplit(names(tab), ","), length)
pick <- names(tab)[sizes >= min_items][1]
if (is.na(pick)) stop("no item block with at least ", min_items, " items")
cols <- as.integer(strsplit(pick, ",")[[1]]); rows <- which(key == pick)
cat("block:", length(cols), "items,", length(rows), "students (",
    round(100 * length(rows) / nrow(cog), 1), "% of the country sample )\n")

# full credit = 1; partial and no credit = 0
Y <- sapply(cols, function(j) as.integer(S[rows, j] == max(S[, j], na.rm = TRUE)))
secs <- Tm[rows, cols]
unit_ms <- median(secs) > 1000                 # item times are stored in milliseconds
if (unit_ms) secs <- secs / 1000
logT <- log(secs)
colnames(Y) <- colnames(logT) <- stems[cols]

q <- read_sav(qqq_file, col_select = c("CNT", "CNTSTUID", "ST004D01T", "ST022Q01TA", "ESCS"))
q <- q[as.character(q$CNT) == cnt, ]
m <- match(as.numeric(cog$CNTSTUID[rows]), as.numeric(q$CNTSTUID))
female <- as.integer(as.numeric(zap_labels(q$ST004D01T[m])) == 1)
ml <- as.integer(as.numeric(zap_labels(q$ST022Q01TA[m])) == 2)
escs <- as.numeric(zap_labels(q$ESCS[m]))

keep <- !is.na(female)                          # G must be observed (it enters the model)
out <- list(Y = Y[keep, ], logT = logT[keep, ], female = female[keep], ml = ml[keep],
            escs = escs[keep], id = as.numeric(cog$CNTSTUID[rows])[keep], items = stems[cols],
            country = cnt, domain = dom,
            log = list(n_country = nrow(cog), n_block = length(rows), n_kept = sum(keep),
                       time_unit = if (unit_ms) "ms" else "s", block = stems[cols]))
f <- file.path(dir, sprintf("prepared_%s_%s.rds", cnt, dom))
saveRDS(out, f)
cat("saved", f, ": n =", sum(keep), ", p =", length(cols),
    "; female", round(mean(out$female), 2), "; ml", round(mean(out$ml, na.rm = TRUE), 2),
    "(NA", sum(is.na(out$ml)), "); escs NA", sum(is.na(out$escs)), "\n")
cat("median item time (s):", round(median(exp(out$logT)), 1),
    "; share of item times < 3 s:", round(mean(exp(out$logT) < 3), 3), "\n")

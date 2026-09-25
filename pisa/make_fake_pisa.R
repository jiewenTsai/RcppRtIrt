# make_fake_pisa.R — synthetic files with the PISA 2018 column layout, to test the pipeline
# without the real data. Two booklet patterns (items 1-10 or 6-15), times in ms, labelled codes;
# multilingual students (ml) are slower on items 1-3 (differential response time).
suppressMessages(library(haven))
dir <- Sys.getenv("PISA_DIR", "data/pisa_fake"); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
set.seed(1); n <- 1600; p <- 15
stems <- sprintf("CM%03dQ01", 1:p)
female <- rbinom(n, 1, .5); ml <- rbinom(n, 1, .25); escs <- rnorm(n) - .4 * ml
th <- rnorm(n) + .3 * escs - .1 * female
tau <- rnorm(n, 0, .4) - .3 * ml
a <- runif(p, .7, 1.5); dd <- rnorm(p, 0, .6); xi <- rnorm(p, 4, .3); g <- rnorm(p, .3, .1)
Y <- matrix(rbinom(n * p, 1, pnorm(outer(th, a) - rep(1, n) %o% dd)), n)
logT <- matrix(xi, n, p, byrow = TRUE) - tau + outer(th, g) + matrix(rnorm(n * p, 0, .5), n)
logT[, 1:3] <- logT[, 1:3] + .4 * ml
book <- rbinom(n, 1, .5)
cog <- data.frame(CNT = "USA", CNTSTUID = 84000000 + 1:n)
for (j in 1:p) {
  adm <- if (j <= 5) book == 0 else if (j >= 11) book == 1 else rep(TRUE, n)
  cog[[paste0(stems[j], "T")]] <- ifelse(adm, round(exp(logT[, j]) * 1000), NA)
  cog[[paste0(stems[j], "S")]] <- labelled_spss(ifelse(adm, Y[, j], 7), c(`Not reached` = 7), na_values = 7)
}
write_sav(cog, file.path(dir, "CY07_MSU_STU_COG.sav"))
qqq <- data.frame(CNT = "USA", CNTSTUID = cog$CNTSTUID,
                  ST004D01T = labelled_spss(ifelse(female == 1, 1, 2), c(Female = 1, Male = 2)),
                  ST022Q01TA = labelled_spss(ifelse(ml == 1, 2, 1), c(`Language of test` = 1, `Other language` = 2)),
                  ESCS = escs)
write_sav(qqq, file.path(dir, "CY07_MSU_STU_QQQ.sav"))
saveRDS(list(theta = th, female = female, ml = ml, escs = escs), file.path(dir, "truth.rds"))
cat("fake PISA files in", dir, "\n")

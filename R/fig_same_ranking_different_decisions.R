# ============================================================
# fig_same_ranking_different_decisions.R  —  SELF-CONTAINED
#
# Core message:
#   Likelihood determines who we think is more likely mosaic (ranking).
#   Prior determines where we draw the line (decision).
#
# Three panels:
#   A. ROC curves under three priors — nearly identical (AUC stable)
#      → ranking is prior-invariant
#
#   B. Posterior distributions for diploid and mosaic embryos
#      under each prior — wildly different location and spread
#      → calibration is prior-dependent
#
#   C. FPR and sensitivity at a single fixed cutoff (0.5)
#      across priors — unstable
#      → decision is prior-dependent
#
# Outputs:
#   figures/fig_same_ranking_different_decisions.png
#   results/prior_dependence_data.csv
# ============================================================

set.seed(2025)
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ── 1. Core functions ─────────────────────────────────────────────────────

SIGMA <- 0.05 * 1.2 / sqrt(10 / 10)   # 0.06

posterior_prob <- function(r_obs, alpha, beta,
                           m_threshold = 0.10, coverage = 10,
                           noise_sd_base = 0.05, pcr_duplicates = 1.2) {
  sigma     <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
  integrand <- function(M)
    dnorm(r_obs, log2(1 + M/2), sigma) * dbeta(M, alpha, beta)
  marg <- tryCatch(
    integrate(integrand, 0, 1, rel.tol=1e-6, abs.tol=1e-9)$value,
    error = function(e) NA_real_)
  if (is.na(marg) || marg <= 0) return(NA_real_)
  numer <- tryCatch(
    integrate(integrand, m_threshold, 1, rel.tol=1e-6, abs.tol=1e-9)$value,
    error = function(e) NA_real_)
  if (is.na(numer)) return(NA_real_)
  min(max(numer / marg, 0), 1)
}

trap_auc <- function(fpr, sens) {
  o <- order(fpr); x <- fpr[o]; y <- sens[o]
  sum(diff(x) * (y[-length(y)] + y[-1]) / 2)
}

# ── 2. Priors ─────────────────────────────────────────────────────────────

PRIORS <- list(
  list(alpha=2,   beta=10, label="Beta(2,10)",
       sublabel="mean\u22480.17\nP(M\u22650.10|prior)=0.70",
       col="#2166AC", lty=1,  pch=19),
  list(alpha=0.5, beta=10, label="Beta(0.5,10)",
       sublabel="mean\u22480.05\nP(M\u22650.10|prior)=0.15",
       col="#1A9850", lty=4,  pch=17),
  list(alpha=1,   beta=49, label="Beta(1,49)",
       sublabel="mean\u22480.02\nP(M\u22650.10|prior)=0.006",
       col="#D73027", lty=5,  pch=15)
)

# ── 3. Simulate data ──────────────────────────────────────────────────────

N      <- 3000
TRUE_M <- 0.15
CUTOFF <- 0.5

cat("Simulating", N, "embryos per class...\n")
r_dip <- rnorm(N, mean = 0,                    sd = SIGMA)
r_mos <- rnorm(N, mean = log2(1 + TRUE_M / 2), sd = SIGMA)

# ── 4. Compute posteriors for all priors ──────────────────────────────────

post_list <- list()

for (pr in PRIORS) {
  cat("Computing posteriors:", pr$label, "\n")
  pd <- vapply(r_dip, posterior_prob, numeric(1),
               alpha=pr$alpha, beta=pr$beta)
  pm <- vapply(r_mos, posterior_prob, numeric(1),
               alpha=pr$alpha, beta=pr$beta)
  post_list[[pr$label]] <- list(diploid=pd, mosaic=pm)
}

# ── 5. Panel A data: ROC curve on likelihood ratio (prior-invariant) ────────
#
# LR(r_obs) = p(r_obs | M = TRUE_M) / p(r_obs | M = 0)
# This is prior-independent: the same LR is computed for all priors.
# Sweeping the LR cutoff gives ONE ROC curve that represents the information
# content of the data — before any prior is applied.
# Overlaying it for each prior (same curve, same colour) makes the point:
# ranking is stable; only the decision boundary moves.

lr_score <- function(r_obs, m = TRUE_M, coverage = 10,
                     noise_sd_base = 0.05, pcr_duplicates = 1.2) {
  sigma <- noise_sd_base * pcr_duplicates / sqrt(coverage / 10)
  # log LR to avoid underflow
  log(dnorm(r_obs, log2(1 + m/2), sigma) + 1e-300) -
  log(dnorm(r_obs, 0,             sigma) + 1e-300)
}

cat("Computing likelihood ratios...
")
lr_dip <- vapply(r_dip, lr_score, numeric(1))
lr_mos <- vapply(r_mos, lr_score, numeric(1))

LR_CUTOFFS <- seq(
  floor(min(c(lr_dip, lr_mos)) * 10) / 10,
  ceiling(max(c(lr_dip, lr_mos)) * 10) / 10,
  by = 0.05
)

nd <- length(lr_dip); nm <- length(lr_mos)
roc_lr <- do.call(rbind, lapply(LR_CUTOFFS, function(co)
  data.frame(cutoff=co,
             fpr  = mean(lr_dip > co),
             sens = mean(lr_mos > co))))
roc_lr$auc <- trap_auc(roc_lr$fpr, roc_lr$sens)
cat(sprintf("  LR-based AUC = %.3f  (same for all priors)
", roc_lr$auc[1]))

# roc_list kept for compatibility but Panel A uses roc_lr
roc_list <- lapply(PRIORS, function(pr) roc_lr)
names(roc_list) <- sapply(PRIORS, `[[`, "label")

# ── 6. Panel C data: FPR and sensitivity at cutoff=0.5 ───────────────────

panel_c <- do.call(rbind, lapply(PRIORS, function(pr) {
  pd  <- post_list[[pr$label]]$diploid
  pm  <- post_list[[pr$label]]$mosaic
  fpr  <- mean(pd > CUTOFF, na.rm=TRUE)
  sens <- mean(pm > CUTOFF, na.rm=TRUE)
  data.frame(prior=pr$label, fpr=fpr, sens=sens,
             prior_mean = pr$alpha/(pr$alpha+pr$beta),
             stringsAsFactors=FALSE)
}))

cat("\nPanel C — at posterior cutoff = 0.5:\n")
print(panel_c)

# ── 7. Save CSV ───────────────────────────────────────────────────────────

write.csv(panel_c, "results/prior_dependence_data.csv", row.names=FALSE)

# ── 8. Plot ───────────────────────────────────────────────────────────────

png("figures/fig_same_ranking_different_decisions.png",
    width=3600, height=1500, res=300)

layout(matrix(1:3, nrow=1), widths=c(1.1, 1.2, 1.0))
par(mar=c(5.5, 5.2, 4.5, 1.0), mgp=c(3.2, 1, 0),
    oma=c(0, 0, 3.0, 0))

# ── Panel A: ROC curves ───────────────────────────────────────────────────

plot(NULL, xlim=c(0,0.3), ylim=c(0,1), las=1,
     xlab="False positive rate",
     ylab="Sensitivity",
     main="A.  ROC curves (posterior cutoff sweep)",
     cex.main=0.92, cex.lab=0.85, cex.axis=0.80)
grid(col="grey93", lty=1)
abline(a=0, b=1, lty=2, col="grey60", lwd=1.0)
text(0.22, 0.19, "random", col="grey50", cex=0.65, srt=50, font=3)

# Single LR-based ROC — prior-invariant
lines(roc_lr$fpr, roc_lr$sens, col="grey20", lwd=2.5, lty=1)

# Mark where each prior's cutoff=0.5 operating point falls on this curve
# by finding the LR cutoff that gives the same FPR as cutoff=0.5 on posterior
for (i in seq_along(PRIORS)) {
  pr   <- PRIORS[[i]]
  ref_fpr <- panel_c$fpr[panel_c$prior == pr$label]
  op   <- roc_lr[which.min(abs(roc_lr$fpr - ref_fpr)), ]
  points(op$fpr, op$sens, pch=pr$pch, cex=1.6, col=pr$col)
  text(op$fpr + 0.012, op$sens - 0.04,
       sprintf("%s
@FPR=%.0f%%
sens=%.0f%%",
               pr$label, op$fpr*100, op$sens*100),
       cex=0.60, col=pr$col, adj=0)
}

legend("bottomright",
       legend=c(sprintf("LR-based ROC  AUC=%.3f", roc_lr$auc[1]),
                sapply(PRIORS, function(p)
                  paste0(p$label, " operating point (posterior cutoff=0.5)"))),
       col=c("grey20", sapply(PRIORS,`[[`,"col")),
       lty=c(1, rep(NA, length(PRIORS))),
       pch=c(NA, sapply(PRIORS,`[[`,"pch")),
       lwd=c(2.5, rep(NA, length(PRIORS))),
       pt.cex=c(NA, rep(1.4, length(PRIORS))),
       bty="n", cex=0.68)

mtext("One ROC curve (likelihood-based)
Three operating points (prior-dependent)",
      side=3, line=0.1, cex=0.68, col="#333333", font=3)

# ── Panel B: posterior distributions ─────────────────────────────────────

par(mar=c(5.5, 4.5, 4.5, 1.0))

# compute density for each prior x class combination
dens_list <- lapply(PRIORS, function(pr) {
  pd <- post_list[[pr$label]]$diploid
  pm <- post_list[[pr$label]]$mosaic
  list(
    dip = density(pd[!is.na(pd)], from=0, to=1, n=512),
    mos = density(pm[!is.na(pm)], from=0, to=1, n=512)
  )
})

# y-axis limit
max_y <- max(sapply(dens_list, function(d)
  max(d$dip$y, d$mos$y))) * 1.05

plot(NULL, xlim=c(0,1), ylim=c(0, max_y), las=1,
     xlab="Posterior P(M \u2265 0.10 | r_obs)",
     ylab="Density",
     main="B.  Posterior distributions per prior",
     cex.main=0.92, cex.lab=0.85, cex.axis=0.80)
grid(col="grey93", lty=1)
abline(v=CUTOFF, lty=2, col="grey30", lwd=1.3)
text(CUTOFF+0.02, max_y*0.97, "cutoff = 0.5",
     cex=0.68, col="grey30", adj=0, font=3)

for (i in seq_along(PRIORS)) {
  pr <- PRIORS[[i]]
  d  <- dens_list[[i]]
  # diploid: solid, lighter
  lines(d$dip$x, d$dip$y,
        col=adjustcolor(pr$col, alpha.f=0.55), lwd=1.6, lty=2)
  # mosaic: solid, full colour
  lines(d$mos$x, d$mos$y,
        col=pr$col, lwd=2.2, lty=pr$lty)
}

legend("topright",
       legend=c(
         sapply(PRIORS, function(p) paste0(p$label, " — mosaic")),
         sapply(PRIORS, function(p) paste0(p$label, " — diploid"))
       ),
       col=c(sapply(PRIORS,`[[`,"col"),
             sapply(PRIORS, function(p) adjustcolor(p$col, 0.55))),
       lty=c(sapply(PRIORS,`[[`,"lty"), rep(2, length(PRIORS))),
       lwd=c(rep(2.2,length(PRIORS)), rep(1.6,length(PRIORS))),
       bty="n", cex=0.62, ncol=2)

mtext("Distributions shift with prior\n\u2192 calibration is prior-dependent",
      side=3, line=0.1, cex=0.68, col="#333333", font=3)

# ── Panel C: bar chart of FPR and sensitivity at cutoff=0.5 ──────────────

par(mar=c(6.5, 5.2, 4.5, 1.5))

n_prior <- nrow(panel_c)
x_pos   <- seq_len(n_prior)
bar_w   <- 0.32
cols    <- sapply(PRIORS, `[[`, "col")

plot(NULL,
     xlim=c(0.5, n_prior+0.5), ylim=c(0, 1),
     xaxt="n", las=1,
     xlab="", ylab="Rate",
     main="C.  FPR and sensitivity at cutoff = 0.5",
     cex.main=0.92, cex.lab=0.85, cex.axis=0.80)
grid(col="grey93", lty=1)
abline(h=c(0.05, 0.10), lty=3, col="grey60")

# FPR bars (left, hatched)
rect(x_pos - bar_w, 0, x_pos, panel_c$fpr,
     col=adjustcolor(cols, alpha.f=0.35),
     border=cols, lwd=1.5, density=25, angle=45)

# Sensitivity bars (right, solid)
rect(x_pos, 0, x_pos + bar_w, panel_c$sens,
     col=adjustcolor(cols, alpha.f=0.75),
     border=cols, lwd=1.5)

# value labels
for (i in x_pos) {
  text(i - bar_w/2, panel_c$fpr[i]  + 0.02,
       sprintf("%.0f%%", panel_c$fpr[i]*100),
       cex=0.72, col=cols[i], font=2)
  text(i + bar_w/2, panel_c$sens[i] + 0.02,
       sprintf("%.0f%%", panel_c$sens[i]*100),
       cex=0.72, col=cols[i], font=2)
}

# x-axis labels with prior info
prior_sublabels <- sapply(PRIORS, function(pr) {
  p_above <- round(1 - pbeta(0.10, pr$alpha, pr$beta), 3)
  paste0(pr$label, "\nP(M\u22650.10|prior)=", p_above)
})
axis(1, at=x_pos + bar_w/2 - bar_w/2,
     labels=prior_sublabels,
     cex.axis=0.68, tick=FALSE, line=0.5)

legend("topright",
       legend=c("FPR (hatched)", "Sensitivity (solid)"),
       fill=c(adjustcolor("grey40",0.35), adjustcolor("grey40",0.75)),
       border="grey40", density=c(25, NA), angle=c(45, NA),
       bty="n", cex=0.75)

mtext("Same cutoff (0.5) \u2192 wildly different FPR\n\u2192 decision is prior-dependent",
      side=3, line=0.1, cex=0.68, col="#333333", font=3)

# ── Overall title ─────────────────────────────────────────────────────────

mtext(
  "Same ranking, different decisions  \u2014  the role of the prior in Bayesian classification",
  outer=TRUE, cex=0.95, font=2, line=1.8)

mtext(
  paste0(
    "Likelihood determines ranking (Panel A: AUC stable across priors).  ",
    "Prior determines calibration (Panel B: posterior distributions shift).  ",
    "A single cutoff (0.5) on posterior produces wildly different FPR ",
    "(Panel C: ", paste(sprintf("%.0f%%", panel_c$fpr*100), collapse=" / "),
    " across priors) \u2014 the decision boundary is prior-dependent."),
  outer=TRUE, cex=0.60, font=3, col="grey25", line=0.3)

dev.off()
cat("\nSaved: figures/fig_same_ranking_different_decisions.png\n")
cat("Saved: results/prior_dependence_data.csv\n")

# ── Console summary ───────────────────────────────────────────────────────

cat("\n=== SUMMARY ===\n\n")
cat(sprintf("Panel A — LR-based AUC = %.3f  (same for all priors)\n",
           roc_lr$auc[1]))

cat(sprintf("\nPanel C — at cutoff = %.1f (decision instability):\n", CUTOFF))
for (i in seq_len(nrow(panel_c)))
  cat(sprintf("  %-15s  FPR = %4.1f%%   sensitivity = %4.1f%%\n",
              panel_c$prior[i],
              panel_c$fpr[i]*100,
              panel_c$sens[i]*100))

cat("\nKey message:\n")
cat(sprintf("  AUC (LR-based, prior-invariant): %.3f\n", roc_lr$auc[1]))
cat("  FPR range:  ",
    sprintf("%.1f%%", min(panel_c$fpr)*100), "—",
    sprintf("%.1f%%", max(panel_c$fpr)*100),
    "  (unstable)\n")
cat("  Sens range: ",
    sprintf("%.1f%%", min(panel_c$sens)*100), "—",
    sprintf("%.1f%%", max(panel_c$sens)*100),
    "  (unstable)\n")

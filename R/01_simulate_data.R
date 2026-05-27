# Synthetic data generator for the salt + CNT stress portfolio piece.
# Two crops (Barley, Cowpea) -- both salt-tolerant-ish C3s with published
# CNT priming literature -- crossed with NaCl (0/50/100/150 mM) and
# multi-walled CNT (0/50 mg/L). 6 replicates per cell.
#
# Effect structure (qualitative, matches the published phenotype space for
# salt-stressed crops):
#   * NaCl: monotonic decrease in growth, monotonic increase in oxidative
#     damage (MDA, H2O2), antioxidant enzymes (SOD, CAT, APX) and proline.
#   * CNT priming: partially rescues growth under salt, lowers damage
#     markers, further upregulates antioxidant enzymes (priming
#     hypothesis).
#   * Cowpea is more salt-sensitive than barley (steeper salt slope on growth).

set.seed(42)

species  <- c("Barley", "Cowpea")
nacl_mM  <- c(0, 50, 100, 150)
cnt_mgL  <- c(0, 50)
n_rep    <- 6

design <- expand.grid(
  species = species,
  nacl    = nacl_mM,
  cnt     = cnt_mgL,
  rep     = seq_len(n_rep),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

# species-specific salt sensitivity (1 = more sensitive)
sens <- c(Barley = 0.7, Cowpea = 1.0)

simulate_row <- function(species, nacl, cnt) {
  s <- sens[[species]]
  # normalised salt dose 0..1
  sd <- nacl / 150
  # CNT rescue factor: 0..1 of the salt-induced loss is recovered
  rescue <- if (cnt > 0) 0.35 else 0
  # CNT also independently boosts enzyme activity (priming)
  prime  <- if (cnt > 0) 0.20 else 0

  # ---- growth ----
  shoot_length <- (28 - 14 * s * sd * (1 - rescue)) *
                  rnorm(1, 1, 0.06)             # cm
  root_length  <- (15 -  9 * s * sd * (1 - rescue)) *
                  rnorm(1, 1, 0.07)
  fresh_weight <- (1.8 - 1.0 * s * sd * (1 - rescue)) *
                  rnorm(1, 1, 0.08)             # g

  # shoot-root axis (hypocotyl + radicle) fresh weight, mg
  shoot_root_axis_weight <- (520 - 280 * s * sd * (1 - rescue)) *
                            rnorm(1, 1, 0.08)   # mg

  # ---- oxidative damage ----
  MDA <- (8  + 22 * s * sd * (1 - 0.6 * rescue)) *
         rnorm(1, 1, 0.10)                       # nmol/g FW
  H2O2 <- (1.2 + 3.0 * s * sd * (1 - 0.6 * rescue)) *
          rnorm(1, 1, 0.10)                      # umol/g FW

  # ---- antioxidant enzymes (salt up, CNT up) ----
  SOD <- (40 + 80 * sd + 25 * prime) * rnorm(1, 1, 0.08)   # U/mg protein
  CAT <- (15 + 28 * sd + 10 * prime) * rnorm(1, 1, 0.09)
  APX <- (12 + 22 * sd +  9 * prime) * rnorm(1, 1, 0.10)

  # ---- osmolyte ----
  proline <- (0.4 + 2.6 * sd + 0.3 * prime) * rnorm(1, 1, 0.12)  # umol/g FW

  list(
    shoot_length_cm = shoot_length,
    root_length_cm  = root_length,
    fresh_weight_g  = fresh_weight,
    shoot_root_axis_weight_mg = shoot_root_axis_weight,
    MDA_nmol_gFW    = MDA,
    H2O2_umol_gFW   = H2O2,
    SOD_U_mg        = SOD,
    CAT_U_mg        = CAT,
    APX_U_mg        = APX,
    proline_umol_gFW= proline
  )
}

vals <- mapply(simulate_row,
               design$species, design$nacl, design$cnt,
               SIMPLIFY = FALSE)
vals <- do.call(rbind.data.frame, lapply(vals, as.data.frame))

out <- cbind(design, vals)
out$treatment <- paste0("NaCl", out$nacl, "_CNT", out$cnt)
out$nacl <- factor(out$nacl, levels = nacl_mM)
out$cnt  <- factor(out$cnt,  levels = cnt_mgL,
                   labels = c("CNT-", "CNT+"))

dir.create("data", showWarnings = FALSE)
write.csv(out, "data/synthetic.csv", row.names = FALSE)
message("Wrote data/synthetic.csv  (", nrow(out), " rows)")

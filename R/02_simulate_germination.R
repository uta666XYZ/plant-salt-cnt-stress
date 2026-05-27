# Synthetic germination time-course for the salt + CNT stress portfolio.
# 2 species (Barley, Cowpea) x 6 NaCl levels (0,50,100,150,200,300 mM)
# x 2 CNT (CNT-/CNT+) x 4 dishes per cell, each dish = 25 seeds.
# Daily counts on Day 0..4. Output: cumulative germination percentage.
#
# Effect structure:
#   * 0..150 mM: nearly complete germination by Day 2-3, slight delay at 150.
#   * 200 mM: partial germination, plateau around 50-60%.
#   * 300 mM: near-zero germination.
#   * CNT+ shifts curves slightly upward / earlier (priming).
#   * Cowpea slightly more salt-sensitive than barley.

set.seed(123)

species <- c("Barley", "Cowpea")
nacl_mM <- c(0, 50, 100, 150, 200, 300)
cnt_lvl <- c("CNT-", "CNT+")
n_dish  <- 4
seeds_per_dish <- 25
days <- 0:4

# species-specific salt sensitivity
sens <- c(Barley = 0.85, Cowpea = 1.0)

# logistic-like germination curve as a function of day
# final asymptote and rate depend on NaCl and CNT
germ_curve <- function(day, asymp, rate, lag) {
  p <- asymp / (1 + exp(-rate * (day - lag)))
  pmax(0, pmin(asymp, p))
}

asymptote <- function(nacl, s, cnt_plus) {
  # baseline asymptote by salt
  base <- dplyr::case_when(
    nacl <=  50 ~ 0.97,
    nacl <= 100 ~ 0.94,
    nacl <= 150 ~ 0.90,
    nacl <= 200 ~ 0.55,
    TRUE        ~ 0.05
  )
  # species sensitivity reduces asymptote at salt
  base <- base - (1 - 1 / s) * 0.05 * (nacl / 300)
  # CNT priming: small lift, more visible at moderate-high salt
  if (cnt_plus) base <- min(0.99, base + 0.04 * (nacl / 300) + 0.01)
  base * (1 - (s - 1) * 0.05 * (nacl > 0))
}

rate_of <- function(nacl, cnt_plus) {
  r <- 3.0 - 1.6 * (nacl / 300)
  if (cnt_plus) r <- r + 0.3
  r
}

lag_of <- function(nacl, cnt_plus) {
  l <- 1.0 + 0.6 * (nacl / 300)
  if (cnt_plus) l <- l - 0.15
  l
}

design <- expand.grid(
  species = species,
  nacl    = nacl_mM,
  cnt     = cnt_lvl,
  dish    = seq_len(n_dish),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

rows <- list()
k <- 0
for (i in seq_len(nrow(design))) {
  sp <- design$species[i]; na <- design$nacl[i]
  cn <- design$cnt[i];     dh <- design$dish[i]
  s  <- sens[[sp]]
  cnt_plus <- cn == "CNT+"
  a <- asymptote(na, s, cnt_plus)
  r <- rate_of(na, cnt_plus)
  L <- lag_of(na, cnt_plus)

  # add small per-dish noise to curve params
  a_d <- max(0, min(1, a + rnorm(1, 0, 0.025)))
  r_d <- max(0.5, r + rnorm(1, 0, 0.15))
  L_d <- max(0.2, L + rnorm(1, 0, 0.1))

  for (d in days) {
    p <- germ_curve(d, a_d, r_d, L_d)
    # binomial draw on 25 seeds, then take CUMULATIVE max so
    # later-day counts cannot drop below earlier-day counts
    if (d == 0) {
      n_germ <- 0
    } else {
      n_germ <- rbinom(1, seeds_per_dish, p)
    }
    k <- k + 1
    rows[[k]] <- data.frame(
      species = sp, nacl = na, cnt = cn, dish = dh,
      day = d, n_germinated = n_germ,
      n_seeds = seeds_per_dish,
      stringsAsFactors = FALSE
    )
  }
}

germ <- do.call(rbind, rows)

# enforce cumulative monotonicity within each (species, nacl, cnt, dish)
suppressMessages({
  if (!requireNamespace("dplyr", quietly = TRUE))
    install.packages("dplyr", repos = "https://cloud.r-project.org")
})
library(dplyr)

germ <- germ %>%
  group_by(species, nacl, cnt, dish) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(n_germinated = cummax(n_germinated),
         germ_pct     = 100 * n_germinated / n_seeds) %>%
  ungroup() %>%
  as.data.frame()

dir.create("data", showWarnings = FALSE)
write.csv(germ, "data/germination.csv", row.names = FALSE)
message("Wrote data/germination.csv  (", nrow(germ), " rows)")

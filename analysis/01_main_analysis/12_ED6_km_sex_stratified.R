# 12_ED6_km_sex_stratified.R — sex-stratified KM, registered as
# ed_sex_stratified_km in panel_registry.R.
#
# Follow-up to 09_F1_km_age.R. Takes the two headline KM analyses:
#   (1) time to any clinical presentation (same cohort as km_presentation)
#   (2) time to rupture, Option 2a strict (same cohort as km_rupture)
# and further stratifies each by sex (Female / Male). Goal is to see
# whether the genotype effect is sex-modifed or whether one sex dominates
# the signal.
#
# For each outcome we produce:
#   • Per-sex KM plot + at-risk table (one ggsurvplot per (sex × endpoint)
#     cell, so each cell ships its own number-at-risk table downstream).
#   • Per-sex log-rank p + per-sex median time-to-event.
#   • Cox PH with genotype + sex + genotype:sex (study + sample_type
#     also adjusted for the rupture model, matching the strict-rupture
#     cohort used by km_rupture in 04b).
#   • Likelihood-ratio test of the genotype × sex interaction term.
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(here)
  library(survival); library(survminer); library(patchwork)
})

source(here("analysis", "helper_scripts", "utils.R"))

df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds")) %>%
  filter(!is.na(mutation_positive)) %>%
  mutate(
    variant_group = case_when(
      !mutation_positive              ~ "Negative",
      mutation == "KRAS G12D"         ~ "KRAS G12D",
      mutation == "KRAS G12V"         ~ "KRAS G12V",
      mutation_gene == "KRAS"         ~ "Other KRAS",
      mutation == "BRAF V600E"        ~ "BRAF V600E",
      mutation_gene == "BRAF"         ~ "Other BRAF",
      TRUE                             ~ "Negative"
    ),
    variant_group = factor(variant_group,
      levels = c("KRAS G12D", "KRAS G12V", "Other KRAS", "BRAF V600E", "Other BRAF", "Negative"))
  )

# PAL_KM (3-arm KM palette) sourced from utils.R via the source() call
# at the top of this script. Restricts to well-powered strata
# (G12D / G12V / Negative) by construction.

# 2026-05-19 (Phase 2 / Iteration 2): ed_sex_stratified_km merged into
# ed_km_diagnostics. The four save_km_panel() calls below now route via
# slot_dir(token) so each panel lands in its canonical
# ed_km_diagnostics/panel_<LETTER>/ slot resolved through the registry.
out_dir   <- here("results", "ExtendedData", "ed_km_diagnostics")
stats_dir <- here("results", "stats")
dir.create(out_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(stats_dir, recursive = TRUE, showWarnings = FALSE)

source(here("analysis", "pipeline", "helpers", "panel_assignments.R"))
.pa_local <- load_panel_assignments()
slot_dir <- function(token) {
  d <- panel_slot_dir(token, .pa_local)
  if (is.null(d)) stop(sprintf(
    "slot_dir: token '%s' is not registered in panel_registry.", token),
    call. = FALSE)
  d
}

source(here("analysis", "pipeline", "helpers", "save_panel.R"))
save_trio <- function(dir, name, plot, w, h)
  save_panel_impl(dir, name, plot, w, h, device = "cairo")

# ─────────────────────────────────────────────────────────────────────────────
# Helper: per-sex KM ggsurvplot (one sex at a time) so we can ship an
# at-risk table with each cell. 2026-05-17: previous implementation used
# ggsurvplot_facet (Female|Male in one ggplot), but that helper returns a
# single ggplot with no $table slot — no at-risk numbers reach Extended
# Data. Switching to two separate ggsurvplot() calls per outcome (Female
# cohort, Male cohort), each with risk.table = TRUE, mirrors the Fig 2
# producer pattern (09_F1_km_age.R) so save_km_panel() can persist curve +
# table separately and compose_figure() can reassemble them downstream.
#
# Returns a ggsurvplot object (NOT a single ggplot); caller passes the
# result straight into save_km_panel().
# ─────────────────────────────────────────────────────────────────────────────
km_sex_panel <- function(dat, time_col, event_col, group_col, sex,
                         xlab, ylab, legend_position = "none") {
  # Use base-R subset to avoid dplyr NSE picking up the wrong `sex`
  # symbol (the function arg vs a potential column reference) — older
  # versions of dplyr have masked the arg here, leaving sub empty.
  sex_val <- sex
  keep <- dat$sex_f == sex_val &
          !is.na(dat[[time_col]]) & !is.na(dat[[event_col]]) &
          !is.na(dat[[group_col]])
  sub <- dat[which(keep), , drop = FALSE]
  sub$gg <- droplevels(factor(sub[[group_col]]))
  # ggsurvplot()'s internals call formula(fit) and expect a fit built via
  # survminer's surv_fit() (which preserves the original formula for
  # ggsurvplot's downstream surv_summary()). Building a formula via
  # sprintf() loses the call attachment, so use `do.call` with a literal
  # formula via bquote so the survfit object carries an inspectable formula.
  form <- bquote(survival::Surv(.(as.name(time_col)), .(as.name(event_col))) ~ gg)
  fit  <- do.call(survminer::surv_fit, list(formula = form, data = sub))

  pal <- as.character(PAL_KM[levels(sub$gg)])

  p <- ggsurvplot(
    legend = legend_position, fit, data = sub,
    # v6.24 (2026-05-20): conf.int.alpha 0.15 -> 0.2 and explicit size
    # = 0.5 added so the sex-stratified curves match Fig 1 KMs D/F/H
    # (linewidth 0.5, CI ribbon alpha 0.2) and the rare_variants KM
    # in this same ED figure 1-to-1.
    conf.int = TRUE, conf.int.alpha = 0.2,
    size      = 0.5,
    palette = pal,
    risk.table = TRUE, risk.table.height = 0.28,
    risk.table.y.text.col = TRUE, risk.table.y.text = FALSE,
    xlab = xlab, ylab = ylab,
    legend.title = "Genotype",
    legend.labs  = display_levels(sub$gg),
    break.time.by = 10, xlim = c(0, KM_AGE_XLIM_MAX),
    censor.size = 1.5,
    ggtheme = theme_avm(),
    tables.theme = theme_cleantable()
  )
  # Dim censor marks — matches 09_F1_km_age.R post-processing.
  p$plot$layers <- lapply(p$plot$layers, function(l) {
    if (inherits(l$geom, "GeomPoint")) l$aes_params$alpha <- 0.4
    l
  })
  # Stamp the sex on each panel so readers don't lose track of which
  # cell is which (panel letters alone don't convey sex).
  p$plot <- p$plot +
    ggtitle(sprintf("%s", sex)) +
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0))
  p
}

# ─────────────────────────────────────────────────────────────────────────────
# OUTCOME 1 — time to any clinical presentation (Fig 4A cohort)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Outcome 1: time to any presentation, sex-stratified ──\n")

km1 <- df %>%
  filter(!is.na(age), !is.na(sex_f),
         variant_group %in% c("KRAS G12D", "KRAS G12V", "Negative")) %>%
  mutate(variant_group = droplevels(variant_group),
         event = 1L, time = age)

cat(sprintf("  n = %d | Female = %d, Male = %d\n",
            nrow(km1), sum(km1$sex_f == "Female"), sum(km1$sex_f == "Male")))
print(km1 %>% count(variant_group, sex_f))

# Per-sex log-rank + per-sex medians
per_sex1 <- lapply(levels(km1$sex_f), function(s) {
  sub <- filter(km1, sex_f == s)
  fit <- survfit(Surv(time, event) ~ variant_group, data = sub)
  lr  <- survdiff(Surv(time, event) ~ variant_group, data = sub)
  list(sex = s, n = nrow(sub),
       p  = 1 - pchisq(lr$chisq, df = length(lr$n) - 1),
       summary = summary(fit)$table)
})
for (r in per_sex1) {
  cat(sprintf("  %s (n=%d): log-rank p = %.3g\n", r$sex, r$n, r$p))
  print(r$summary)
}

# Cox with genotype × sex interaction (no study/sample_type for presentation —
# every patient presents, no censoring, so covariate set matches Fig 4A
# framing which is unadjusted).
cox1_main <- coxph(Surv(time, event) ~ variant_group + sex_f,            data = km1)
cox1_int  <- coxph(Surv(time, event) ~ variant_group * sex_f,            data = km1)
lrt1      <- anova(cox1_main, cox1_int, test = "Chisq")
cat("\n  Cox main-effects model:\n");     print(summary(cox1_main)$coefficients)
cat("\n  Cox interaction model:\n");      print(summary(cox1_int)$coefficients)
cat(sprintf("\n  LRT genotype × sex: Chisq = %.3f, df = %d, p = %.3g\n",
            lrt1$Chisq[2], lrt1$Df[2], lrt1$`Pr(>|Chi|)`[2]))

# Per-sex panels with at-risk tables. Female panel keeps the legend
# (top of figure); Male panel suppresses it so the bar of strata labels
# doesn't duplicate.
# Suppress per-panel legends on every cell so the long "Panel-negative"
# label doesn't get clipped in the narrower 2x2 composite. Composer 46
# adds a single figure-wide Genotype legend strip below the grid.
p_pres_F <- km_sex_panel(km1, "time", "event", "variant_group",
                         sex = "Female",
                         xlab = "Age (years)",
                         ylab = "Proportion symptom-free",
                         legend_position = "none")
p_pres_M <- km_sex_panel(km1, "time", "event", "variant_group",
                         sex = "Male",
                         xlab = "Age (years)",
                         ylab = "Proportion symptom-free",
                         legend_position = "none")
save_km_panel(slot_dir("km_pres_sex_F"), "km_pres_sex_F", p_pres_F, w = 10, h = 8)
save_km_panel(slot_dir("km_pres_sex_M"), "km_pres_sex_M", p_pres_M, w = 10, h = 8)

# ─────────────────────────────────────────────────────────────────────────────
# OUTCOME 2 — time to rupture (Option 2a strict — same cohort as km_rupture)
# ─────────────────────────────────────────────────────────────────────────────
cat("\n── Outcome 2: time to rupture, sex-stratified ──\n")

km2 <- df %>%
  mutate(km2_event = case_when(
    rupture_category == "Ruptured at presentation" ~ 1L,
    rupture_category == "Never ruptured"            ~ 0L,
    TRUE                                            ~ NA_integer_)) %>%
  filter(!is.na(km2_event), !is.na(age), !is.na(sex_f),
         variant_group %in% c("KRAS G12D", "KRAS G12V", "Negative")) %>%
  mutate(variant_group = droplevels(variant_group))

cat(sprintf("  n = %d | events = %d\n",
            nrow(km2), sum(km2$km2_event == 1)))
print(km2 %>% count(variant_group, sex_f, km2_event))

per_sex2 <- lapply(levels(km2$sex_f), function(s) {
  sub <- filter(km2, sex_f == s)
  fit <- survfit(Surv(age, km2_event) ~ variant_group, data = sub)
  lr  <- survdiff(Surv(age, km2_event) ~ variant_group, data = sub)
  list(sex = s, n = nrow(sub),
       events = sum(sub$km2_event == 1),
       p  = 1 - pchisq(lr$chisq, df = length(lr$n) - 1),
       summary = summary(fit)$table)
})
for (r in per_sex2) {
  cat(sprintf("  %s (n=%d, events=%d): log-rank p = %.3g\n",
              r$sex, r$n, r$events, r$p))
  print(r$summary)
}

# Adjusted Cox (matches km_rupture covariates: sex + study + sample_type),
# with and without genotype × sex interaction.
cox2_df <- km2 %>%
  filter(!is.na(study_clean), !is.na(sample_type_clean))
cox2_main <- coxph(Surv(age, km2_event) ~ variant_group + sex_f +
                     study_clean + sample_type_clean, data = cox2_df)
cox2_int  <- coxph(Surv(age, km2_event) ~ variant_group * sex_f +
                     study_clean + sample_type_clean, data = cox2_df)
lrt2 <- anova(cox2_main, cox2_int, test = "Chisq")
cat(sprintf("\n  Cox n = %d\n", nrow(cox2_df)))
cat("\n  Cox main-effects model:\n");    print(summary(cox2_main)$coefficients)
cat("\n  Cox interaction model:\n");     print(summary(cox2_int)$coefficients)
cat(sprintf("\n  LRT genotype × sex: Chisq = %.3f, df = %d, p = %.3g\n",
            lrt2$Chisq[2], lrt2$Df[2], lrt2$`Pr(>|Chi|)`[2]))

p_rupt_F <- km_sex_panel(km2, "age", "km2_event", "variant_group",
                         sex = "Female",
                         xlab = "Age (years)",
                         ylab = "Proportion rupture-free",
                         legend_position = "none")
p_rupt_M <- km_sex_panel(km2, "age", "km2_event", "variant_group",
                         sex = "Male",
                         xlab = "Age (years)",
                         ylab = "Proportion rupture-free",
                         legend_position = "none")
save_km_panel(slot_dir("km_rupt_sex_F"), "km_rupt_sex_F", p_rupt_F, w = 10, h = 8)
save_km_panel(slot_dir("km_rupt_sex_M"), "km_rupt_sex_M", p_rupt_M, w = 10, h = 8)

# ─────────────────────────────────────────────────────────────────────────────
# Write stats dump
# ─────────────────────────────────────────────────────────────────────────────
fmt_per_sex <- function(lst, ev = FALSE) {
  unlist(lapply(lst, function(r) {
    hd <- if (ev)
      sprintf("### %s (n=%d, events=%d) — log-rank p = %.3g",
              r$sex, r$n, r$events, r$p)
    else
      sprintf("### %s (n=%d) — log-rank p = %.3g", r$sex, r$n, r$p)
    c(hd, capture.output(print(r$summary)), "")
  }))
}

writeLines(c(
  sprintf("# %s — KM analyses further stratified by sex",
          panel_prose_tag(c("km_pres_sex_F", "km_pres_sex_M",
                            "km_rupt_sex_F", "km_rupt_sex_M"))),
  sprintf("Generated %s", Sys.Date()),
  "",
  sprintf("## %s — Time to any clinical presentation",
          panel_prose_tag(c("km_pres_sex_F", "km_pres_sex_M"))),
  sprintf("Cohort n = %d (Female %d, Male %d); everyone presents (event=1).",
          nrow(km1), sum(km1$sex_f == "Female"), sum(km1$sex_f == "Male")),
  "",
  fmt_per_sex(per_sex1, ev = FALSE),
  "### Cox PH — main effects (genotype + sex)",
  capture.output(print(summary(cox1_main)$coefficients)),
  "",
  "### Cox PH — with genotype × sex interaction",
  capture.output(print(summary(cox1_int)$coefficients)),
  "",
  sprintf("LRT genotype × sex: Chisq = %.3f, df = %d, p = %.3g",
          lrt1$Chisq[2], lrt1$Df[2], lrt1$`Pr(>|Chi|)`[2]),
  "",
  sprintf("## %s — Time to rupture (Option 2a strict)",
          panel_prose_tag(c("km_rupt_sex_F", "km_rupt_sex_M"))),
  sprintf("Cohort n = %d (events = %d).",
          nrow(km2), sum(km2$km2_event == 1)),
  "",
  fmt_per_sex(per_sex2, ev = TRUE),
  "### Cox PH — main effects (genotype + sex + study + sample_type)",
  capture.output(print(summary(cox2_main)$coefficients)),
  "",
  "### Cox PH — with genotype × sex interaction",
  capture.output(print(summary(cox2_int)$coefficients)),
  "",
  sprintf("LRT genotype × sex: Chisq = %.3f, df = %d, p = %.3g",
          lrt2$Chisq[2], lrt2$Df[2], lrt2$`Pr(>|Chi|)`[2])
), file.path(stats_dir, "km_sex_stratified_stats.txt"))

# ─────────────────────────────────────────────────────────────────────────────
# Stats manifest fragment — edfig22 sex-stratified KM
# -----------------------------------------------------------------------------
# §4 prose cites "genotype × sex interaction P = 0.712" for the presentation
# outcome (cox1 LRT). The rupture outcome LRT is also emitted as the prose
# sensitivity note in ED22 caption ("directionally consistent but formally
# underpowered at 3 df").
# ─────────────────────────────────────────────────────────────────────────────
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))

.per_sex_get <- function(per_sex_list, sex, field) {
  for (r in per_sex_list) {
    if (r$sex == sex) return(r[[field]])
  }
  NA
}

edfig22_fragment <- list(
  # Outcome 1 (time to presentation)
  n_pres               = nrow(km1),
  n_pres_female        = sum(km1$sex_f == "Female"),
  n_pres_male          = sum(km1$sex_f == "Male"),
  logrank_p_pres_female = .per_sex_get(per_sex1, "Female", "p"),
  logrank_p_pres_male   = .per_sex_get(per_sex1, "Male",   "p"),
  cox_lrt_chisq_pres   = unname(lrt1$Chisq[2]),
  cox_lrt_df_pres      = unname(lrt1$Df[2]),
  cox_lrt_p_pres       = unname(lrt1$`Pr(>|Chi|)`[2]),
  # Outcome 2 (time to rupture)
  n_rupt               = nrow(km2),
  n_rupt_events        = sum(km2$km2_event == 1L),
  n_rupt_female        = .per_sex_get(per_sex2, "Female", "n"),
  n_rupt_male          = .per_sex_get(per_sex2, "Male",   "n"),
  events_rupt_female   = .per_sex_get(per_sex2, "Female", "events"),
  events_rupt_male     = .per_sex_get(per_sex2, "Male",   "events"),
  logrank_p_rupt_female = .per_sex_get(per_sex2, "Female", "p"),
  logrank_p_rupt_male   = .per_sex_get(per_sex2, "Male",   "p"),
  cox_lrt_chisq_rupt   = unname(lrt2$Chisq[2]),
  cox_lrt_df_rupt      = unname(lrt2$Df[2]),
  cox_lrt_p_rupt       = unname(lrt2$`Pr(>|Chi|)`[2]),
  # Full Cox coefficient tables
  cox_main_coef_pres   = as.data.frame(summary(cox1_main)$coefficients),
  cox_int_coef_pres    = as.data.frame(summary(cox1_int)$coefficients),
  cox_main_coef_rupt   = as.data.frame(summary(cox2_main)$coefficients),
  cox_int_coef_rupt    = as.data.frame(summary(cox2_int)$coefficients)
)
write_stats_section(section = "edfig22", stats = edfig22_fragment)

cat("\nWrote:\n  ", file.path(out_dir, "km_pres_sex_{F,M}.{pdf,png,rds}"),
    "\n  ", file.path(out_dir, "km_rupt_sex_{F,M}.{pdf,png,rds}"),
    "\n  ", file.path(stats_dir, "km_sex_stratified_stats.txt"), "\n")

# 24_F4_dpcr_waterfall.R — Fig 4 panel d: endovascular liquid-biopsy dPCR.
#
# Real QIAcuity nanoplate (26K) partition-amplitude "waterfall" for the GREEN
# (FAM / KRAS G12D) channel, four sample roles in story order:
#   KRAS G12D+ control  (FFPE tumour positive ctrl, Run 3 well C1)
#   Peri-nidal cfDNA    (endovascular liquid biopsy,        well F1)
#   Peripheral cfDNA    (systemic comparator,               well G1)
#   Germline            (buffy-coat control,                well H2)
#
# Input: data/raw/UAB Priority Samples Run 3 …/…_G_… (this repo's raw QIAcuity
#   GREEN-channel export, read directly — provenance fix 2026-05-27; the biopsy
#   data must not be sourced from the public single-cell sibling repo).
#   Gate = 100 RFU; the instrument Is-positive call == RFU > gate.
#   NTC (H3) is the negative control for ED Fig 3, not this panel.
#   Schema: role, role_label, partition_id, fluorescence_rfu, threshold, call.
#
# Output: results/Figure4/biopsy_dpcr_waterfall/biopsy_dpcr_waterfall.{png,pdf,rds}
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(ggplot2); library(here); library(ggtext)
})
source(here("analysis", "helper_scripts", "utils.R"))
source(here("analysis", "pipeline", "helpers", "save_panel.R"))

out_root <- here("results", "Figure4")
save_panel <- function(subdir, name, plot, w, h)
  save_panel_impl(file.path(out_root, subdir), name, plot, w, h, device = "cairo")

# Provenance: read the internal raw QIAcuity export from THIS repo's data/raw/
# (the biopsy data must not be sourced from the public single-cell sibling repo).
# Fig 4d is the KRAS G12D (FAM / "G") channel across four wells of Run 3.
raw_dir <- here("data", "raw",
  "UAB Priority Samples Run 3 011626_RFU_img1_G_Y_O_R_C_Fr_26_05_2026_09_32_12_UTC-04_00")
g_csv <- list.files(raw_dir, pattern = "img1_G_[0-9].*\\.csv$", full.names = TRUE)
if (length(g_csv) != 1L)
  stop("Fig 4d raw KRAS G12D channel export not found in ", raw_dir)

role_levels <- c("KRAS G12D+ control", "Peri-nidal cfDNA",
                 "Peripheral cfDNA", "Germline")
# Well -> role for this run: C1 G12D+ control, F1 peri-nidal cfDNA (AVMUAB019),
# G1 peripheral cfDNA (AVMUAB039), H2 germline (AVMUAB019). Story order = role_levels.
fig4_wells <- c(C1 = "KRAS G12D+ control", F1 = "Peri-nidal cfDNA",
                G1 = "Peripheral cfDNA",  H2 = "Germline")
# ggtext-markdown facet labels (italic gene + superscript variant for the ctrl).
role_labeller <- ggplot2::as_labeller(c(
  "KRAS G12D+ control" = "*KRAS*<sup>G12D+</sup> control",
  "Peri-nidal cfDNA"   = "Peri-nidal cfDNA",
  "Peripheral cfDNA"   = "Peripheral cfDNA",
  "Germline"           = "Germline"))

# `skip = 1` drops the QIAcuity "sep=," preamble line; positivity is the
# instrument call at the applied threshold; invalid partitions (Is positive = NA)
# are excluded, matching the prior extract.
dat <- read_csv(g_csv, skip = 1, show_col_types = FALSE) %>%
  filter(Well %in% names(fig4_wells), !is.na(`Is positive`)) %>%
  transmute(role_label = factor(unname(fig4_wells[Well]), levels = role_levels),
            fluorescence_rfu = RFU,
            threshold = Threshold,
            call = factor(ifelse(`Is positive` == 1, "positive", "negative"),
                          levels = c("negative", "positive"))) %>%
  group_by(role_label) %>% mutate(part_x = row_number()) %>% ungroup()

gate <- unique(dat$threshold)
stopifnot(length(gate) == 1L)

# Downsample negatives for a light vector PDF; keep every positive.
set.seed(42)
neg_keep <- dat %>% filter(call == "negative") %>%
  group_by(role_label) %>% slice_sample(n = 2000) %>% ungroup()
plot_df <- bind_rows(dat %>% filter(call == "positive"), neg_keep) %>%
  arrange(call)   # negatives first, positives on top

# Positive partitions are coloured by the dye channel (FAM = green emission),
# matching ED Fig 3; see PAL_FLUOROPHORE in utils.R.
mut_col <- PAL_FLUOROPHORE[["FAM"]]
pal <- c(negative = "#D0D0D0", positive = mut_col)
sz  <- c(negative = 0.30,      positive = 0.65)
al  <- c(negative = 0.40,      positive = 0.90)

p <- ggplot(plot_df, aes(part_x, fluorescence_rfu, colour = call,
                         size = call, alpha = call)) +
  geom_point(stroke = 0) +
  ref_hline(yintercept = gate, kind = "threshold") +
  facet_wrap(~ role_label, nrow = 1, labeller = role_labeller) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_size_manual(values = sz, guide = "none") +
  scale_alpha_manual(values = al, guide = "none") +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.01, 0.10))) +
  scale_x_continuous(breaks = c(0, 10000, 20000),
                     labels = c("0", "10,000", "20,000"),
                     expand = expansion(mult = c(0.03, 0.03))) +
  labs(x = "Partition", y = "Fluorescence (RFU)") +
  theme_nature_panel() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position  = "none",
        strip.text = element_markdown(size = NM$body_pt, face = "bold",
                                      family = NM$font_family))

save_panel("biopsy_dpcr_waterfall", "biopsy_dpcr_waterfall", p, w = 5.0, h = 2.1)
message("[fig4_dpcr_waterfall] wrote panel to ",
        file.path(out_root, "biopsy_dpcr_waterfall"))

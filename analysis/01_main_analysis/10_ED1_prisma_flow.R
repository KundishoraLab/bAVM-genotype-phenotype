# 10_ED1_prisma_flow.R
# ─────────────────────────────────────────────────────────────────────────────
# Renders the PRISMA 2020 literature-screening waterfall for the published-
# arm systematic review (registry token `ed_prisma`, Extended Data Fig. 1).
#
# Counts are final (verified 2026-05-23):
#   Identification:  1,948 records via PubMed/MEDLINE
#                    4 removed before screening (2 duplicates + 2 automation)
#   Screening:       1,944 screened → 1,928 excluded → 16 sought
#                    0 not retrieved → 16 assessed for eligibility
#   Eligibility:     11 excluded (wrong design n=9, wrong pop n=1,
#                    wrong sample size n=1)
#   Included:        5 studies → 186 patient-level bAVM genotype-phenotype pairs
#
# Output: results/ExtendedData/ed_prisma/prisma_flow.{png,pdf,svg}
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(ggplot2); library(here); library(tibble); library(dplyr); library(grid)
})

source(here("analysis", "helper_scripts", "utils.R"))

# ── Final counts ─────────────────────────────────────────────────────────────
N_IDENTIFIED  <- 1948
N_DEDUP       <- 2
N_AUTOMATION  <- 2
N_SCREENED    <- 1944   # = 1948 - 2 - 2
N_EXCL_SCREEN <- 1928
N_SOUGHT      <- 16
N_NOT_RETR    <- 0
N_ASSESSED    <- 16
N_EXCL_DESIGN <- 9
N_EXCL_POP    <- 1
N_EXCL_SIZE   <- 1
N_INCLUDED    <- 5
N_PATIENTS    <- 186

# ── Box geometry ─────────────────────────────────────────────────────────────
# Main column at x = 0.  Exclusion column at x = 3.6 (no overlap).
# Arrows branch rightward off the vertical waterfall at the midpoint
# between successive main boxes (PRISMA 2020 convention).
#
#   Main box:  half-width 1.5 → right edge at x = 1.5
#   Excl box:  center 3.6, half-width 1.4 → left edge at x = 2.2
#   H-arrow:   from x = 1.5 → x = 2.2  (0.7-unit gap, arrow points right)
MW  <- 1.5   # main-box half-width
MH  <- 0.50  # main-box half-height
EX  <- 3.6   # exclusion column centre
EW  <- 1.4   # excl-box half-width

main_boxes <- tribble(
  ~id,           ~y,    ~label,
  "identified",   7.2,  sprintf("Records identified from\nPubMed/MEDLINE\n(n = %s)",
                                 format(N_IDENTIFIED, big.mark = ",")),
  "screened",     5.4,  sprintf("Records screened\n(n = %s)",
                                 format(N_SCREENED, big.mark = ",")),
  "sought",       3.6,  sprintf("Reports sought for retrieval\n(n = %d)", N_SOUGHT),
  "assessed",     1.8,  sprintf("Reports assessed for eligibility\n(n = %d)", N_ASSESSED),
  "included",     0.0,  sprintf("Studies included in review\n(n = %d)\n%d patient-level bAVM\ngenotype–phenotype pairs",
                                 N_INCLUDED, N_PATIENTS)
) %>% mutate(x = 0, w = MW, h = MH)

excl_boxes <- tribble(
  ~id,             ~y,    ~eh,   ~label,
  "removed",        6.30,  0.60,  sprintf("Records removed before screening:\nDuplicate records (n = %d)\nMarked ineligible by automation (n = %d)",
                                           N_DEDUP, N_AUTOMATION),
  "excl_screen",    4.50,  0.40,  sprintf("Records excluded\n(n = %s)",
                                           format(N_EXCL_SCREEN, big.mark = ",")),
  "not_retr",       2.70,  0.35,  sprintf("Reports not retrieved\n(n = %d)", N_NOT_RETR),
  "excl_elig",      0.90,  0.60,  sprintf("Reports excluded:\nWrong study design (n = %d)\nWrong population (n = %d)\nWrong sample size (n = %d)",
                                           N_EXCL_DESIGN, N_EXCL_POP, N_EXCL_SIZE)
) %>% mutate(x = EX, w = EW)

# ── Arrows ───────────────────────────────────────────────────────────────────
# Vertical waterfall: tip at top edge of lower box.
v_arrows <- tibble(
  x_from = 0, y_from = main_boxes$y[-nrow(main_boxes)] - MH,
  x_to   = 0, y_to   = main_boxes$y[-1] + MH
)

# Horizontal branches: FROM right edge of main column TO left edge of excl box.
# Arrow points rightward → arrowhead is on the excl-box side.
h_arrows <- tibble(
  x_from = MW,             # right edge of main box
  y_from = excl_boxes$y,
  x_to   = EX - EW,       # left edge of exclusion box
  y_to   = excl_boxes$y
)

# ── Stage labels (left margin) ────────────────────────────────────────────────
stage_labels <- tribble(
  ~y,    ~label,
  7.2,   "Identification",
  4.5,   "Screening",
  1.8,   "Eligibility",
  0.0,   "Included"
)

# ── Build plot ───────────────────────────────────────────────────────────────
p <- ggplot() +
  # Main waterfall boxes
  geom_rect(data = main_boxes,
            aes(xmin = x - w, xmax = x + w,
                ymin = y - h,  ymax = y + h),
            fill = "white", colour = "grey30", linewidth = 0.5) +
  geom_text(data = main_boxes,
            aes(x = x, y = y, label = label),
            size = NM$text$body_mm, lineheight = 1.1, colour = "#1a1a1a") +
  # Exclusion boxes
  geom_rect(data = excl_boxes,
            aes(xmin = x - w, xmax = x + w,
                ymin = y - eh, ymax = y + eh),
            fill = "white", colour = "grey30", linewidth = 0.5) +
  geom_text(data = excl_boxes,
            aes(x = x, y = y, label = label),
            size = NM$text$body_mm, lineheight = 1.1, colour = "#1a1a1a") +
  # Vertical waterfall arrows
  geom_segment(data = v_arrows,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               colour = "grey35", linewidth = 0.4) +
  # Horizontal branch arrows
  geom_segment(data = h_arrows,
               aes(x = x_from, y = y_from, xend = x_to, yend = y_to),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               colour = "grey35", linewidth = 0.4) +
  # Stage labels in left margin
  geom_text(data = stage_labels,
            aes(x = -2.2, y = y, label = label),
            hjust = 0, fontface = "bold",
            size = NM$text$body_mm, colour = "grey20") +
  coord_cartesian(xlim = c(-2.4, 5.2), ylim = c(-0.8, 7.9), clip = "off") +
  theme_void() +
  theme(plot.margin = margin(12, 12, 12, 40))

# ── Save ─────────────────────────────────────────────────────────────────────
out_dir <- here("results", "ExtendedData", "ed_prisma")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "prisma_flow.pdf"), p,
       width = 7.20, height = 6.00, device = cairo_pdf, family = NM$font_family)
ggsave(file.path(out_dir, "prisma_flow.png"), p,
       width = 7.20, height = 6.00, dpi = 300, type = "cairo")
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(out_dir, "prisma_flow.svg"), p,
         width = 7.20, height = 6.00, device = svglite::svglite,
         fix_text_size = FALSE)
}

readme <- file.path(out_dir, "README.md")
if (file.exists(readme)) file.remove(readme)

cat(sprintf("✓ PRISMA literature-screening flow written: %s\n",
            file.path(out_dir, "prisma_flow.png")))

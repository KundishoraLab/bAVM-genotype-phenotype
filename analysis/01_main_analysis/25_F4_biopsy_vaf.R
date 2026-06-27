# 25_F4_biopsy_vaf.R — TERT-referenced VAF for the Fig 4d endovascular
# liquid-biopsy case, computed from THIS repo's raw QIAcuity export.
# -----------------------------------------------------------------------------
# Provenance: the biopsy dPCR data lives here in data/raw/ (it must NOT be
# sourced from the public single-cell repo). This producer reads the raw
# per-partition channel exports directly.
#
# Method (Andy's, and consistent with the FFPE cohort VAFs): VAF = mutant-KRAS
# positive partitions / TERT positive partitions, at the QIAcuity-applied
# thresholds recorded in the export. TERT (Cy5.5) is the assay's DNA-input
# reference (see Methods); KRAS (12p) and TERT (5p) are single-copy autosomal
# loci, so under copy-number neutrality TERT alleles ≈ total KRAS alleles,
# making mutant/TERT an estimate of VAF when no wild-type-KRAS probe exists.
#
# Fig 4d wells (this run): C1 = G12D+ control, F1 = peri-nidal cfDNA
# (AVMUAB019), G1 = peripheral cfDNA (AVMUAB039), H2 = germline (AVMUAB019).
# The matched resected-lesion VAF is the cohort tissue value for AVMUAB019,
# computed by the same TERT-referenced method (already in the master VAF
# column); cfDNA and tissue are independent measurements (genotype-concordant,
# VAFs not expected to be identical).
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(here)
})
source(here("analysis", "pipeline", "helpers", "write_stats_section.R"))

raw_dir <- here("data", "raw",
  "UAB Priority Samples Run 3 011626_RFU_img1_G_Y_O_R_C_Fr_26_05_2026_09_32_12_UTC-04_00")

read_channel <- function(ch) {
  f <- list.files(raw_dir, pattern = sprintf("img1_%s_[0-9].*\\.csv$", ch),
                  full.names = TRUE)
  if (length(f) != 1L) stop(sprintf("expected 1 '%s' channel CSV, found %d", ch, length(f)))
  d <- fread(f, skip = 1, sep = ",", header = TRUE, showProgress = FALSE)
  setnames(d, make.names(names(d)))
  d
}

PERINIDAL_WELL <- "F1"   # AVMUAB019 peri-nidal cfDNA (see well map above)

g  <- read_channel("G")    # KRAS G12D (FAM)
fr <- read_channel("Fr")   # TERT (Cy5.5), DNA-input reference

g_pos <- g[Well == PERINIDAL_WELL,  sum(Is.positive == 1, na.rm = TRUE)]
t_pos <- fr[Well == PERINIDAL_WELL, sum(Is.positive == 1, na.rm = TRUE)]
# VAF = mutant / (mutant + TERT) copies, matching the QIAcuity Mutation
# Detection formula documented in Methods (TERT as the wild-type/total-allele
# reference). mutant << TERT, so this ≈ mutant/TERT.
cfdna_vaf <- 100 * g_pos / (g_pos + t_pos)

# Matched resected-lesion (tissue) VAF — cohort value for the same patient.
df <- readRDS(here("data", "processed", "bAVM_analysis_ready.rds"))
tissue_row <- df %>% filter(grepl("AVMUAB019", patient_id, ignore.case = TRUE))
if (nrow(tissue_row) != 1L)
  warning(sprintf("expected 1 AVMUAB019 tissue row, found %d", nrow(tissue_row)))

write_stats_section(section = "biopsy", stats = list(
  cfdna_patient         = "AVMUAB019",
  cfdna_g12d_partitions = as.integer(g_pos),
  cfdna_tert_partitions = as.integer(t_pos),
  cfdna_vaf_pct         = round(cfdna_vaf, 4),
  tissue_sample_uid     = tissue_row$sample_uid[1],
  tissue_vaf_pct        = tissue_row$vaf_pct[1]
))

cat(sprintf("── Fig 4d biopsy VAF ──\n  cfDNA (F1, AVMUAB019): G12D %d / TERT %d = %.2f%%\n  tissue (%s): %.2f%%\n",
            g_pos, t_pos, cfdna_vaf, tissue_row$sample_uid[1], tissue_row$vaf_pct[1]))

# =============================================================================
# compose_figure.R — native-size submission-composite builder
# -----------------------------------------------------------------------------
# Producer scripts save per-panel ggplot OBJECTS as RDS alongside their
# standalone PDF/PNG. compose_figure() reads those RDS objects, applies a
# uniform native-typography theme overlay, lays them out via patchwork at the
# manuscript footprint specified in NATIVE (utils.R), and ggsaves at exactly
# that footprint with NO publisher-stage resize. The result is a single
# composite PNG + PDF that already matches its target column width and
# height in print, so axis text lands inside Nature Medicine's 5–7pt spec
# directly off the producer rather than relying on the journal scaling
# 14pt text down.
#
# Inputs:
#   panel_rds    Named list of paths to per-panel RDS files. Names are panel
#                tokens (registry tokens), values are absolute paths. Each
#                RDS must contain a ggplot or patchwork object.
#   layout       patchwork design string (areas given by capital letters).
#                Letter -> token mapping is provided via `area_to_token`.
#   area_to_token  Named character vector mapping single-letter area codes
#                in `layout` to panel tokens (which in turn key into
#                `panel_rds`). Allows the layout shape to be independent of
#                the manifest's logical panel order.
#   out_stem     Output path stem (no extension). Composer writes
#                <stem>.pdf via cairo_pdf and <stem>.png via cairo at 300dpi.
#   width_in,
#   height_in    Native render dimensions in inches. Width should match a
#                Nature Medicine column target (NATIVE$widths_in).
#   base_size    ggplot base_size to overlay on every panel (defaults to
#                NATIVE$base_size). Set per-figure if a particular composite
#                needs slightly different scaling.
#   panel_tag_size  pt size for the bold a/b/c panel letter; defaults to
#                NATIVE$panel_tag.
#
# Returns invisibly the assembled patchwork object (for inspection in tests).
#
# Conventions:
#   * Every panel RDS contains a ready-to-render ggplot. compose_figure()
#     does NOT re-fit data, only re-themes for typography uniformity.
#   * If a panel is itself a patchwork (e.g. ggsurvplot's curve+table block),
#     pass the result of patchwork::wrap_elements(full = ...) at producer
#     time so the composer treats it as a single grob.
#   * panel_letter() in caption fns reads from the registry-driven cache;
#     compose_figure adds a matching `tag_levels` annotation so the rendered
#     letter and the caption agree by construction.
#
# Composite dimensions (rationale for choosing width_in / height_in):
#
#   Nature Medicine prints at 7.2 in double-column width. Producers render
#   above that and the publisher scales down at typeset time, which lets
#   us keep base_size 14 panel typography on disk and still land within
#   the 5-7 pt journal text floor.
#
#   - **14 in width** is the project default (Fig 1, Fig 3, ED4-12, etc.).
#     Scales 2x down to 7.2 in print. Each cell of a 2-col grid resolves
#     to ~7 in standalone -> ~3.5 in print. KM curves with at-risk tables,
#     forests, and dumbbells all render legibly at this density.
#   - **16 in width** (Fig 2, ED14): used when 3-col grids need each cell
#     at ~5 in standalone (~2.5 in print) without packing too tightly.
#     Scales 2.2x down.
#   - **22 in width** (Fig 4 only): the EC-subtype + variant + size legend
#     strip on the right margin (guides = "collect") consumes ~4 in;
#     bumping the composite to 22 in keeps each of the 3 panel columns at
#     a full ~6 in even after the legend strip lands. Scales 3.05x down.
#
#   Per-row heights follow native panel aspect: KM panels need ~5-6 in
#   for the curve + at-risk table, dumbbells / waffles read fine at ~4 in,
#   and stacked bars read fine at ~3 in. ROW_HEIGHTS in each composer
#   should be tuned to the tallest panel in each row. Don't shrink to fit
#   the figure on a screen at the cost of crowding the print rendering.
#
# Idempotent — safe to re-run.
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(patchwork)
  # ggtext is loaded so element_grob.element_markdown is dispatched correctly
  # for any panel whose axis labels were authored with ggtext::element_markdown
  # (italic gene names, super/subscripts). Without this, the composer's
  # ggsave call dies with "no applicable method for 'element_grob' applied to
  # an object of class 'element_markdown'".
  library(ggtext)
})

source(here::here("analysis", "helper_scripts", "utils.R"))

#' Pass-through panel adapter for the composer.
#'
#' compose_figure() is a pure-layout helper: it does NOT re-theme panels.
#' Each producer is responsible for the typography of its own panels — we
#' just lay them out in a grid. The only transformation we make here is
#' wrapping patchwork sub-composites (e.g. ggsurvplot curve + risk-table
#' stacks) via `wrap_elements(full = ...)` so the parent patchwork sees
#' them as a single cell instead of trying to flatten their internal grobs
#' into separate design-string slots.
.apply_native_theme <- function(p, base_size = NATIVE$base_size) {
  if (inherits(p, "patchwork")) return(patchwork::wrap_elements(full = p))
  p
}

#' Read a panel input file and adapt it for the parent patchwork.
#'
#' Four input shapes are supported:
#'   1. `<panel>.png` — static raster (e.g. CONSORT flowchart from
#'      DiagrammeR / external schematics). Read via png::readPNG +
#'      grid::rasterGrob, wrapped as a single cell.
#'   2. `<panel>.rds` containing
#'      `list(kind = "km_panel", curve, table, ratio)` — written by
#'      `save_km_panel()`. Curve and table reassembled at the saved ratio
#'      and wrapped as one cell. Producers own typography.
#'   3. `<panel>.rds` containing a patchwork — wrapped as a single grob.
#'   4. `<panel>.rds` containing a ggplot — passed through unchanged.
#'
#' @param path absolute path to the producer-written panel file.
#' @param post hook applied to the raw ggplot before any wrap_elements()
#'   coercion happens. Used by compose_figure() to inject
#'   `+ theme(legend.position = "none")` when strip_legends = TRUE, which
#'   must run on the raw ggplot so the theme propagates — once a panel
#'   becomes a wrap_elements() grob its theme is frozen.
#' @return ggplot or patchwork object ready for layout.
.read_panel_rds <- function(path, base_size = NATIVE$base_size, post = identity) {
  if (!file.exists(path)) {
    stop(sprintf("[compose_figure] panel file missing: %s\n  (producer script must save the panel as RDS or PNG)", path),
         call. = FALSE)
  }

  # Static raster panel (e.g. CONSORT flowchart). Wrap as grob so patchwork
  # treats it as a single cell — the figure-level composer doesn't need to
  # know the producer pipeline (DiagrammeR, Illustrator, etc.).
  if (tolower(tools::file_ext(path)) == "png") {
    img <- png::readPNG(path)
    g   <- grid::rasterGrob(img, interpolate = TRUE)
    return(patchwork::wrap_elements(full = g))
  }

  obj <- readRDS(path)

  # KM-panel metadata pointer (from save_km_panel) — rebuild from
  # components at the producer's saved height ratio. Producers own the
  # typography of curve and table; we just reassemble the block here so
  # patchwork treats it as a single cell. `post` is applied to the curve
  # ggplot (not the table) so legend-stripping affects the curve while
  # leaving the at-risk table's typography intact.
  if (is.list(obj) && identical(obj[["kind"]], "km_panel")) {
    panel_dir <- dirname(path)
    curve <- post(readRDS(file.path(panel_dir, obj$curve)))
    tab   <- readRDS(file.path(panel_dir, obj$table))
    block <- curve / tab + patchwork::plot_layout(heights = obj$ratio)
    return(patchwork::wrap_elements(full = block))
  }

  .apply_native_theme(post(obj), base_size = base_size)
}

compose_figure <- function(panel_rds,
                           layout,
                           area_to_token,
                           out_stem,
                           width_in,
                           height_in,
                           base_size      = NATIVE$base_size,
                           panel_tag_size = NATIVE$panel_tag,
                           heights        = NULL,
                           widths         = NULL,
                           guides         = NULL,
                           strip_legends  = FALSE,
                           bottom_strip   = NULL,
                           bottom_strip_height = 0.6,
                           side_strip     = NULL,
                           side_strip_width = 1.8,
                           native_theme   = FALSE,
                           outer_margin_in = 0.3,
                           dpi            = 300) {
  # `strip_legends`: when TRUE, applies `& theme(legend.position = "none")` to
  # the composite, suppressing every per-panel legend AFTER the design is
  # laid out. Use when panels carry differently-shaped legends that would
  # otherwise consume varying amounts of cell space and break row/column
  # plot-region alignment (canonical case: Fig 2's 3×3 grid where some
  # panels show a Genotype legend and others show a Score legend). Pair
  # with `bottom_strip` to add a single figure-level legend bar that
  # documents the palettes.
  #
  # `bottom_strip`: optional ggplot/patchwork rendered as an additional row
  # below the main grid. Typically a horizontally-arranged set of mini
  # legend plots (see 26_F1_assemble.R for a worked example). Height in
  # inches is controlled by `bottom_strip_height`; that height is added on
  # top of `height_in` so the main grid keeps its allocated geometry.
  #
  # `side_strip`: optional ggplot/patchwork attached as an additional COLUMN
  # to the right of the main grid (and below the bottom_strip, if any).
  # Width in inches is controlled by `side_strip_width`; that width is added
  # to `width_in` so the main grid keeps its allocated horizontal geometry.
  # Use for a vertically-stacked legend column when the figure wants its
  # palette guides off to one side rather than below. Both bottom_strip
  # and side_strip can be supplied simultaneously — bottom is appended
  # first (vertically), then side is appended next to that composite.
  #
  # Locking panel sizes: pass `heights` and `widths` as `grid::unit(...,
  # "in")` to force absolute panel dimensions. With numeric (relative)
  # values, patchwork distributes the available space proportionally. With
  # absolute units, each design-grid cell renders at exactly the specified
  # inch size regardless of figure-level decorations — the recommended
  # mode going forward so panels stay locked across composites.
  #
  # `native_theme`: when TRUE, applies theme_avm_native() (utils.R, base 9
  # for Nature Medicine print size) to every loaded panel via the `post`
  # hook. Producers save panels at theme_avm() base 18 for legible
  # standalone PNGs at native ggsave size; composer overlays base 9 so the
  # composite renders correctly when the composite footprint is at
  # journal-column dimensions. Without this, 18pt text inside a 4.5-in
  # composite cell looks oversized and clips against cell boundaries.
  # `guides`: forwarded to patchwork::plot_layout(guides = ...). Project
  # convention:
  #
  #   "collect" — use when ≥2 panels share a palette (same scale +
  #               same factor levels). Patchwork merges the duplicate
  #               legends and pulls them to a single side strip on
  #               the figure margin, freeing horizontal panel area.
  #               Examples: Fig 3 (dumbbells C/D share BINARY_PALETTE),
  #               Fig 4 (UMAPs C + stacked bar E + lollipop F all use
  #               EC_SUBTYPE_COLORS), ED10 (6 panels share Gene legend),
  #               ED15 (per-patient UMAPs share atlas palette).
  #
  #   NULL (default) — use when each panel carries a unique aesthetic
  #               or panels suppress their own legends explicitly.
  #               Patchwork keeps each panel's legend in place.
  #               Examples: ED11 per-variant phenotype (panels B/C
  #               suppress legends; A keeps its own forest legend),
  #               figures whose panels have visually distinct palettes.
  #
  # Mixed cases (one panel needs a unique inline legend that should NOT
  # collect) are typically handled by setting `show.legend = FALSE` on
  # the per-panel guide that should stay inline, or by composing that
  # panel as a sub-patchwork with `guides = "keep"` before passing it
  # in. See 33_ED9_assemble.R for a worked example.
  stopifnot(is.list(panel_rds), !is.null(names(panel_rds)),
            is.character(layout), length(layout) == 1L,
            is.character(area_to_token), !is.null(names(area_to_token)),
            is.character(out_stem), length(out_stem) == 1L)

  # Validate every area letter in the layout has a token mapping AND the
  # token resolves to a panel RDS. Fail fast with a useful message.
  area_codes <- unique(unlist(strsplit(gsub("[^A-Z#]", "", layout), "")))
  area_codes <- setdiff(area_codes, "#")  # `#` is patchwork's empty-cell sentinel
  missing_areas <- setdiff(area_codes, names(area_to_token))
  if (length(missing_areas) > 0L) {
    stop(sprintf("[compose_figure] layout uses areas with no token mapping: %s",
                 paste(missing_areas, collapse = ", ")), call. = FALSE)
  }
  needed_tokens <- unique(unname(area_to_token[area_codes]))
  missing_tokens <- setdiff(needed_tokens, names(panel_rds))
  if (length(missing_tokens) > 0L) {
    stop(sprintf("[compose_figure] layout needs panels with no RDS path: %s",
                 paste(missing_tokens, collapse = ", ")), call. = FALSE)
  }

  # Read panels in area-letter order so the patchwork additions and the
  # layout's design string agree (patchwork associates the i-th `+`-added
  # plot with the i-th area letter when sorted alphabetically).
  # `post` hook applies legend-stripping to the RAW ggplot before
  # .read_panel_rds() potentially wraps the panel in wrap_elements(), which
  # would freeze its theme and prevent legend.position = "none" from taking
  # effect. (Current patchwork class hierarchy makes every ggplot also
  # inherit "patchwork", so the inherits() check inside .apply_native_theme
  # always wraps — meaning post-wrap theme injection is a no-op in practice.)
  # Compose the post hook from `strip_legends` and `native_theme` flags.
  # Both use `&` (patchwork's "apply to all children") rather than `+`.
  # Saved panels are stored as S7-patchwork objects (every ggplot now
  # subclasses patchwork), and `+ theme(...)` on such an object adds the
  # theme at the OUTER patchwork level only — the inner ggplot keeps its
  # original theme settings. The `&` operator descends through patchwork
  # children and modifies each child plot directly.
  .apply_native <- isTRUE(native_theme)
  .strip_legs  <- isTRUE(strip_legends)
  post <- function(p) {
    if (.apply_native) p <- p & theme_avm_native(base_size = base_size)
    if (.strip_legs)  p <- p & theme(legend.position = "none")
    # Re-assert plot.tag after theme overrides. theme_classic (called
    # inside theme_avm_native) resets plot.tag to NULL when applied via
    # patchwork's `&`, which kills the panel-letter tag set by
    # plot_annotation(tag_levels = "a") below. Re-applying it here keeps
    # tags visible regardless of whether native_theme is on.
    p <- p & theme(
      plot.tag = element_text(face = "bold", size = panel_tag_size)
    )
    p
  }
  area_codes_sorted <- sort(area_codes)
  panels <- lapply(area_codes_sorted, function(code) {
    tok <- area_to_token[[code]]
    val <- panel_rds[[tok]]
    # Inline panels — when the panel_rds value is already a ggplot/patchwork
    # object (built in the composer rather than loaded from a producer RDS
    # file), use it directly. Skip the post hook so author-built legend
    # columns / annotation panels are not legend-stripped by strip_legends.
    if (!is.character(val)) return(val)
    .read_panel_rds(val, base_size = base_size, post = post)
  })

  # Build the patchwork composite. Reduce starts with the first panel and
  # `+`s the rest in sorted area order, then applies plot_layout(design = ...)
  # so the design string maps each area letter to the correctly-ordered slot.
  comp <- Reduce(`+`, panels[-1L], init = panels[[1L]])
  comp <- comp +
    plot_layout(design = layout,
                heights = heights,
                widths  = widths,
                guides  = guides) +
    plot_annotation(
      tag_levels = "a",
      theme = theme(
        plot.tag = element_text(face = "bold", size = panel_tag_size),
        # plot.margin sets the figure-level breathing room around the
        # entire composite. `outer_margin_in` (inches) is converted to pt
        # here AND added to the ggsave canvas size below — without that
        # canvas bump, the margin steals space from the design grid and
        # patchwork compresses cells to fit, clipping axis labels.
        plot.margin = margin(outer_margin_in, outer_margin_in,
                             outer_margin_in, outer_margin_in, unit = "in")
      )
    )

  # Optional alignment helpers — applied AFTER the design layout so they
  # affect every cell of the main grid uniformly. The bottom_strip is
  # stacked OUTSIDE the design (as a second-level patchwork) so it does
  # not participate in the grid's row/column unification and therefore
  # cannot push panel regions around.
  total_height <- height_in
  total_width  <- width_in
  # `strip_legends` is handled in the per-panel `post` hook above; no
  # post-assembly `&` operator is needed (and wouldn't work anyway, since
  # most panels are already wrap_elements grobs by the time the composite
  # is assembled).
  if (!is.null(bottom_strip)) {
    # CRITICAL: the main composite is wrapped via wrap_elements(full = comp)
    # BEFORE the `/` join. Without that wrap, attaching a bottom_strip
    # discards the inner plot_annotation(tag_levels = "a") — patchwork
    # rebuilds the outer composite and the tag info attached to the inner
    # composite is lost. wrap_elements freezes the inner composite as an
    # opaque grob whose tags are already baked into its rendered grobs,
    # so they survive the outer layout. bottom_strip is also wrapped so
    # it's not tagged "b" by the outer composite's tag system.
    comp <- patchwork::wrap_elements(full = comp) /
            patchwork::wrap_elements(full = bottom_strip) +
      plot_layout(heights = grid::unit(c(height_in, bottom_strip_height), "in"))
    total_height <- height_in + bottom_strip_height
  }
  if (!is.null(side_strip)) {
    # Same wrap-then-join pattern as bottom_strip — preserves inner tags.
    # Parens are load-bearing: `+` binds tighter than `|`, so without them
    # `plot_layout(...)` would attach to wrap_elements(side_strip) only
    # rather than to the combined left/right composite.
    comp <- (patchwork::wrap_elements(full = comp) |
             patchwork::wrap_elements(full = side_strip)) +
      plot_layout(widths = grid::unit(c(total_width, side_strip_width), "in"))
    total_width <- total_width + side_strip_width
  }

  # ggsave at exact native dimensions. cairo backend so unicode glyphs (β,
  # ≤, ≥, ρ) render correctly in the embedded fonts.
  #
  # Set the same seed before each ggsave so geom_jitter() (and any other
  # random-position layers) lay their dots in the SAME positions in the
  # PDF and PNG renders. Without this, the two outputs diverge for any
  # panel whose content includes jittered points (Fig 1C VAF violins,
  # ED VAF×phenotype, etc.) — collaborators flipping between formats see
  # the same data at different jitter offsets, which is visually
  # confusing. Hash of `out_stem` keeps the seed stable per composite
  # but distinct across composites.
  .seed <- abs(sum(utf8ToInt(basename(out_stem)))) %% .Machine$integer.max
  dir.create(dirname(out_stem), recursive = TRUE, showWarnings = FALSE)
  set.seed(.seed)
  # Canvas dimensions = design-grid totals + outer margin on every side.
  # Without the margin bump, plot.margin (which sits OUTSIDE the grid)
  # would steal space from the cells, compressing them and clipping axis
  # labels at the figure edges.
  canvas_w <- total_width  + 2 * outer_margin_in
  canvas_h <- total_height + 2 * outer_margin_in
  ggsave(paste0(out_stem, ".pdf"), comp,
         width = canvas_w, height = canvas_h,
         device = cairo_pdf, family = NM$font_family)
  set.seed(.seed)
  ggsave(paste0(out_stem, ".png"), comp,
         width = canvas_w, height = canvas_h,
         dpi = dpi, type = "cairo")
  # v6.19 (2026-05-20): emit SVG alongside PDF/PNG so every composite is
  # round-trippable into Adobe Illustrator with native text. svglite is
  # listed in DESCRIPTION; failure to load is soft-warned, not fatal.
  if (requireNamespace("svglite", quietly = TRUE)) {
    set.seed(.seed)
    tryCatch(
      ggsave(paste0(out_stem, ".svg"), comp,
             width = canvas_w, height = canvas_h,
             device = svglite::svglite, fix_text_size = FALSE),
      error = function(e) warning(sprintf(
        "compose_figure: SVG export failed for %s — %s",
        out_stem, conditionMessage(e)), call. = FALSE)
    )
  }

  message(sprintf("  ✓ composite: %s.{pdf,png,svg} (%.2f × %.2f in @ base_size=%d)",
                  out_stem, canvas_w, canvas_h, base_size))
  invisible(comp)
}

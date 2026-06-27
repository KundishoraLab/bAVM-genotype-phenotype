# =============================================================================
# resolve_panels.R — prose-side panel/table token resolver
# -----------------------------------------------------------------------------
# Consumes analysis/pipeline/panel_registry.R + narrative .qmd files containing
# tokens of the form:
#     @fig[tok]            -> **Fig. Nx**
#     @fig[tok1,tok2]      -> **Fig. Na,b**       (2 panels, comma)
#     @fig[tok1,tok2,tok3] -> **Fig. Na-c**       (3+ contiguous, hyphen)
#     @fig[tok1,tok3]      -> **Fig. Na,c**       (non-contiguous, comma)
#     @edfig[...]          -> **Extended Data Fig. ...**  (same rules)
#     @tab[tok]            -> **Table N**
#     @supptab[tok1,tok2]  -> **Supplementary Tables N, M**
#
# Derived numbering rules:
#   - Figure numbers per track: assigned in order a group's FIRST referenced
#     token appears in prose (walking narrative sections in main_text.qmd
#     include order).
#   - Panel letters per group: assigned in order a panel token FIRST appears.
#     This means composite layout should follow prose order (or the composite
#     builder should honour the resolver-assigned letters).
#   - Table numbers per track: assigned in order a table token FIRST appears.
#
# Contiguity collapsing inside one macro call:
#   1 token                     -> "x"
#   2 tokens                    -> "x, y"  (always comma, never hyphen)
#   >=3 consecutive tokens      -> "x-z"   (hyphen range)
#   Mixed                       -> "a-c, e"  (ranges where possible)
#
# Entry points (intended use):
#   source("analysis/pipeline/helpers/resolve_panels.R")
#   resolution <- resolve_sections(section_paths, registry_path, sections_root)
#   resolution$text[[path]]   # resolved string
#   resolution$maps           # token/group/table assignment diagnostics
#
# Phase 1 scope: returns resolved text in-memory. Writing resolved .qmd files
# or hooking into run_all.R is intentionally deferred.
# =============================================================================

# ---- registry loading ------------------------------------------------------

load_panel_registry <- function(path = NULL) {
  if (is.null(path)) {
    path <- file.path(here::here(), "analysis", "pipeline", "panel_registry.R")
  }
  # Source into a fresh env so we do not leak globals. Parent is baseenv
  # so assignment and list() resolve; user code still cannot see our symbols.
  env <- new.env(parent = baseenv())
  sys.source(path, envir = env)
  if (!exists("panel_registry", envir = env, inherits = FALSE)) {
    stop("panel_registry not defined in ", path, call. = FALSE)
  }
  env$panel_registry
}

# Build token -> (track, group, group_panels) index and validate uniqueness.
build_token_index <- function(registry) {
  idx <- list()

  for (grp_name in names(registry$figures)) {
    grp <- registry$figures[[grp_name]]
    if (is.null(grp$track) || is.null(grp$panels)) {
      stop(sprintf("figures[[%s]] must have $track and $panels", grp_name),
           call. = FALSE)
    }
    if (!grp$track %in% c("Fig", "EDFig")) {
      stop(sprintf("figures[[%s]]$track must be 'Fig' or 'EDFig' (got '%s')",
                   grp_name, grp$track), call. = FALSE)
    }

    # Register the group NAME as its own citation handle. `@fig[group_name]`
    # resolves to "**Fig. N**" (no panel letter) — useful when prose cites
    # a whole figure rather than specific panels.
    if (grp_name %in% names(idx)) {
      stop(sprintf("group name '%s' collides with an existing token",
                   grp_name), call. = FALSE)
    }
    idx[[grp_name]] <- list(
      kind         = "group",
      track        = grp$track,
      group        = grp_name,
      group_panels = grp$panels
    )

    # Register each panel token.
    for (tok in grp$panels) {
      if (tok %in% names(idx)) {
        stop(sprintf("duplicate panel token '%s' (seen in %s and %s)",
                     tok, idx[[tok]]$group, grp_name), call. = FALSE)
      }
      idx[[tok]] <- list(
        kind         = "panel",
        track        = grp$track,
        group        = grp_name,
        group_panels = grp$panels
      )
    }
  }

  for (tok in names(registry$tables)) {
    track <- registry$tables[[tok]]
    if (!track %in% c("Table", "SuppTable")) {
      stop(sprintf("tables[[%s]] must be 'Table' or 'SuppTable' (got '%s')",
                   tok, track), call. = FALSE)
    }
    if (tok %in% names(idx)) {
      stop(sprintf("token '%s' collides between figures and tables", tok),
           call. = FALSE)
    }
    idx[[tok]] <- list(kind = "table", track = track, group = NULL)
  }

  idx
}

# ---- text loader -----------------------------------------------------------

.load_text <- function(path) {
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

# ---- citation scan ---------------------------------------------------------

# Returns a data-frame-like list with one row per macro call, each row carrying
# its section index, character offset, macro name, and token list.
parse_citations <- function(txt, section_ord) {
  pat <- "@(fig|edfig|tab|supptab)\\[([^\\]]+)\\]"
  m   <- gregexpr(pat, txt, perl = TRUE)[[1]]
  if (length(m) == 1L && m[1] == -1L) {
    return(list(rows = list(), text = txt))
  }
  hits <- regmatches(txt, gregexpr(pat, txt, perl = TRUE))[[1]]
  rows <- vector("list", length(hits))
  for (i in seq_along(hits)) {
    cg <- regmatches(hits[i], regexec(pat, hits[i], perl = TRUE))[[1]]
    macro   <- cg[2]
    raw     <- cg[3]
    tokens  <- trimws(strsplit(raw, ",", fixed = TRUE)[[1]])
    tokens  <- tokens[nzchar(tokens)]
    rows[[i]] <- list(
      macro       = macro,
      tokens      = tokens,
      section_ord = section_ord,
      pos         = m[i],
      length      = attr(m, "match.length")[i],
      match_text  = hits[i]
    )
  }
  list(rows = rows, text = txt)
}

# ---- macro -> track mapping ------------------------------------------------

.macro_to_track <- function(macro) {
  switch(macro,
    fig      = "Fig",
    edfig    = "EDFig",
    tab      = "Table",
    supptab  = "SuppTable",
    stop("unknown macro: ", macro, call. = FALSE)
  )
}

.track_label_inline <- function(track) {
  switch(track,
    Fig       = "Fig.",
    EDFig     = "Extended Data Fig.",
    Table     = "Table",
    SuppTable = "Supplementary Table",
    track
  )
}

.track_label_plural <- function(track) {
  switch(track,
    SuppTable = "Supplementary Tables",
    Table     = "Tables",
    track  # Figs plural rarely needed inline; fall through as-is
  )
}

# ---- range collapsing ------------------------------------------------------

# Collapse a sorted integer vector into a "1-3, 5, 7-9"-style string.
collapse_ints <- function(ints) {
  ints <- sort(unique(ints))
  if (!length(ints)) return("")
  if (length(ints) == 1L) return(as.character(ints[1]))
  # Split into runs of consecutive integers.
  breaks <- c(1L, which(diff(ints) != 1L) + 1L)
  runs <- split(ints, cumsum(seq_along(ints) %in% breaks))
  parts <- vapply(runs, function(r) {
    if (length(r) >= 3L) sprintf("%d-%d", r[1], r[length(r)])
    else paste(r, collapse = ", ")
  }, character(1L))
  paste(parts, collapse = ", ")
}

collapse_letters <- function(letter_vec) {
  # Dispatch through integer positions, then map back to letters.
  ords <- sort(unique(match(letter_vec, letters)))
  if (!length(ords)) return("")
  if (length(ords) == 1L) return(letters[ords[1]])
  breaks <- c(1L, which(diff(ords) != 1L) + 1L)
  runs <- split(ords, cumsum(seq_along(ords) %in% breaks))
  parts <- vapply(runs, function(r) {
    if (length(r) >= 3L) sprintf("%s-%s", letters[r[1]], letters[r[length(r)]])
    else paste(letters[r], collapse = ", ")
  }, character(1L))
  paste(parts, collapse = ", ")
}

# ---- assignment pass -------------------------------------------------------

# Walk citations in document order and assign:
#   - group_number[track][group]  <- integer
#   - panel_letter[group][token]  <- character (a, b, c, ...)
#   - table_number[track][token]  <- integer
# Errors on unknown tokens or cross-group mixing in one macro call.
assign_numbers <- function(citations, token_index) {
  group_number <- list(Fig = integer(0L), EDFig = integer(0L))
  panel_letter <- list()
  table_number <- list(Table = integer(0L), SuppTable = integer(0L))

  for (c in citations) {
    track_from_macro <- .macro_to_track(c$macro)

    # Resolve each token.
    if (length(c$tokens) == 0L) {
      stop(sprintf("empty token list in @%s[] at section %d pos %d",
                   c$macro, c$section_ord, c$pos), call. = FALSE)
    }
    infos <- lapply(c$tokens, function(tok) {
      if (!tok %in% names(token_index)) {
        stop(sprintf("unknown token '%s' in @%s[...] at section %d pos %d",
                     tok, c$macro, c$section_ord, c$pos), call. = FALSE)
      }
      token_index[[tok]]
    })

    # Validate track consistency with macro.
    for (i in seq_along(infos)) {
      if (infos[[i]]$track != track_from_macro) {
        stop(sprintf("token '%s' is %s but cited via @%s[...] (expects %s)",
                     c$tokens[i], infos[[i]]$track, c$macro, track_from_macro),
             call. = FALSE)
      }
    }

    if (track_from_macro %in% c("Fig", "EDFig")) {
      # All tokens in one @fig[...]/@edfig[...] must share a group.
      groups <- unique(vapply(infos, `[[`, character(1L), "group"))
      if (length(groups) != 1L) {
        stop(sprintf("@%s[%s] mixes panels from groups (%s); cite separately",
                     c$macro, paste(c$tokens, collapse = ","),
                     paste(groups, collapse = ", ")), call. = FALSE)
      }
      grp <- groups[1]

      # Disallow mixing a group-name reference with individual panel tokens
      # inside one macro call — the author should pick one style.
      kinds <- unique(vapply(infos, `[[`, character(1L), "kind"))
      if (length(kinds) != 1L) {
        stop(sprintf(
          "@%s[%s] mixes a group reference with panel tokens; use one or the other",
          c$macro, paste(c$tokens, collapse = ",")), call. = FALSE)
      }

      # Assign group number on first sighting in this track (applies to
      # both group-name refs and panel-token refs).
      if (!grp %in% names(group_number[[track_from_macro]])) {
        group_number[[track_from_macro]][[grp]] <-
          length(group_number[[track_from_macro]]) + 1L
      }

      # Panel-letter assignment only applies when citing panels (not when
      # citing the figure as a whole).
      if (kinds == "panel") {
        if (is.null(panel_letter[[grp]])) {
          panel_letter[[grp]] <- character(0L)
        }
        for (tok in c$tokens) {
          if (!tok %in% names(panel_letter[[grp]])) {
            next_idx <- length(panel_letter[[grp]]) + 1L
            if (next_idx > 26L) {
              stop(sprintf("group '%s' exceeded 26 panels", grp), call. = FALSE)
            }
            panel_letter[[grp]][[tok]] <- letters[next_idx]
          }
        }
      }

    } else {
      # Table track.
      for (tok in c$tokens) {
        if (!tok %in% names(table_number[[track_from_macro]])) {
          table_number[[track_from_macro]][[tok]] <-
            length(table_number[[track_from_macro]]) + 1L
        }
      }
    }
  }

  # --- Fill-in pass: uncited panels in an active group ----------------------
  # A group is "active" once it has received a figure NUMBER, either via a
  # panel citation (assigns number + letter) or via a group-name citation
  # like @edfig[<group_name>] (assigns number only). For every active group,
  # ensure all registry-listed panels have letter assignments by walking the
  # registry in order. Letters already set by prose citation keep their
  # prose-order assignment; remaining panels get the next letters in
  # registry order.
  #
  # This keeps caption functions (which reference every panel by token) and
  # composite builders (which need a label per slot) working even when the
  # prose cites only the figure as a whole (e.g. "Extended Data Fig. 5"
  # rather than "Extended Data Fig. 5a,b,c,d").
  #
  # First-citation ordering for the check is unaffected: the check only
  # parses actual prose citations, not these derived letters.
  active_groups <- unique(c(names(group_number[["Fig"]]),
                            names(group_number[["EDFig"]])))
  for (grp in active_groups) {
    if (is.null(panel_letter[[grp]])) panel_letter[[grp]] <- character(0L)
    grp_tokens_in_registry <- character(0L)
    for (tok in names(token_index)) {
      info <- token_index[[tok]]
      if (!is.null(info$group) && identical(info$group, grp) &&
          identical(info$kind, "panel")) {
        grp_tokens_in_registry <- c(grp_tokens_in_registry, tok)
      }
    }
    uncited <- setdiff(grp_tokens_in_registry, names(panel_letter[[grp]]))
    for (tok in uncited) {
      next_idx <- length(panel_letter[[grp]]) + 1L
      if (next_idx > 26L) break
      panel_letter[[grp]][[tok]] <- letters[next_idx]
    }
  }

  list(group_number = group_number,
       panel_letter = panel_letter,
       table_number = table_number)
}

`%||%` <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.atomic(a) && length(a) == 1L && is.na(a)) return(b)
  a
}

# ---- rendering one citation -------------------------------------------------

render_one <- function(c, maps) {
  track <- .macro_to_track(c$macro)

  if (track %in% c("Fig", "EDFig")) {
    # Look up the first token to determine group + citation kind
    # (group-name vs panel).
    info <- .token_index_cache[[c$tokens[1]]]
    grp  <- info$group
    n    <- maps$group_number[[track]][[grp]]

    if (info$kind == "group") {
      # Whole-figure citation, no panel letter.
      body <- sprintf("%s %d", .track_label_inline(track), n)
    } else {
      # Panel citation(s). Collapse contiguous letters per the range rule;
      # single-panel groups and full-figure citations (all panels cited
      # together) render with no letter.
      lets <- vapply(c$tokens, function(tok) maps$panel_letter[[grp]][[tok]],
                     character(1L))
      all_group_panels <- get_group_panels(grp)
      all_cited <- length(c$tokens) == length(all_group_panels) &&
                   setequal(c$tokens, all_group_panels)
      if ((length(lets) == 1L && is_single_panel_group(grp)) || all_cited) {
        body <- sprintf("%s %d", .track_label_inline(track), n)
      } else {
        body <- sprintf("%s %d%s", .track_label_inline(track), n,
                        collapse_letters(lets))
      }
    }
  } else {
    # Tables: number(s) only.
    nums <- vapply(c$tokens, function(tok) maps$table_number[[track]][[tok]],
                   integer(1L))
    if (length(nums) == 1L) {
      body <- sprintf("%s %d", .track_label_inline(track), nums[1])
    } else {
      body <- sprintf("%s %s", .track_label_plural(track), collapse_ints(nums))
    }
  }
  sprintf("**%s**", body)
}

# Helpers that need access to registry panel lists; we stash a reference into
# maps during assignment.
find_group_for_token <- function(tok, maps) {
  for (grp in names(maps$panel_letter)) {
    if (tok %in% names(maps$panel_letter[[grp]])) return(grp)
  }
  NULL
}
get_group_panels <- function(grp) {
  # Panel token list lives in the registry; we pull it at resolve time.
  panels <- .registry_cache$figures[[grp]]$panels
  if (is.null(panels)) character(0L) else panels
}
is_single_panel_group <- function(grp) {
  length(get_group_panels(grp)) == 1L
}

# Module-level caches so render_one can see the registry + token index
# without threading them through every argument. Set by resolve_sections().
.registry_cache    <- NULL
.token_index_cache <- NULL

# ---- main entry point ------------------------------------------------------

# section_paths: character vector of full paths, in rendered-document order.
# registry: optional pre-loaded list; loaded from default path if NULL.
resolve_sections <- function(section_paths,
                             registry = NULL,
                             registry_path = NULL) {
  if (is.null(registry)) registry <- load_panel_registry(registry_path)
  .registry_cache <<- registry
  on.exit(.registry_cache <<- NULL, add = TRUE)

  token_index <- build_token_index(registry)
  .token_index_cache <<- token_index
  on.exit(.token_index_cache <<- NULL, add = TRUE)

  # --- Pass 1: parse all citations in document order ---
  all_citations <- list()
  sec_text <- character(length(section_paths))
  names(sec_text) <- section_paths
  for (i in seq_along(section_paths)) {
    txt <- .load_text(section_paths[i])
    sec_text[i] <- txt
    parsed <- parse_citations(txt, section_ord = i)
    all_citations <- c(all_citations, parsed$rows)
  }

  # --- Pass 2: sort by document order (section then pos) and assign ---
  if (length(all_citations) > 0L) {
    ord <- order(
      vapply(all_citations, `[[`, integer(1L), "section_ord"),
      vapply(all_citations, `[[`, numeric(1L), "pos")
    )
    all_citations <- all_citations[ord]
  }
  maps <- assign_numbers(all_citations, token_index)

  # Citation-order gate — warn (not stop) on ascending-order violations.
  check_citation_order(all_citations, maps, token_index)

  # --- Pass 3: rewrite each section's text ---
  resolved <- sec_text
  for (i in seq_along(section_paths)) {
    txt <- sec_text[i]
    if (!grepl("@(fig|edfig|tab|supptab)\\[", txt, perl = TRUE)) next
    # Re-parse this section alone so we have correct offsets (they are
    # section-local in parse_citations).
    parsed <- parse_citations(txt, section_ord = i)
    if (length(parsed$rows) == 0L) next
    # Walk right-to-left so character offsets remain valid.
    ord_desc <- order(vapply(parsed$rows, `[[`, numeric(1L), "pos"),
                      decreasing = TRUE)
    for (j in ord_desc) {
      c <- parsed$rows[[j]]
      replacement <- render_one(c, maps)
      start <- c$pos
      end   <- c$pos + c$length - 1L
      txt <- paste0(substr(txt, 1L, start - 1L),
                    replacement,
                    substr(txt, end + 1L, nchar(txt)))
    }
    resolved[i] <- collapse_adjacent_fig_refs(txt)
  }

  list(
    text        = resolved,
    maps        = maps,
    token_index = token_index
  )
}

# ---- Post-processing: collapse adjacent whole-figure refs -------------------
# Runs of adjacent same-track whole-figure citations separated by "; " are
# collapsed using the canonical range rules (comma for 2, hyphen for 3+
# consecutive). Only fires on refs WITHOUT panel letters — e.g.
#   **Extended Data Fig. 1**; **Extended Data Fig. 2**  ->  **Extended Data Figs. 1, 2**
#   **Extended Data Fig. 1**; **Extended Data Fig. 2**; **Extended Data Fig. 3**  ->  **Extended Data Figs. 1-3**
# "Extended Data Fig." is processed before "Fig." to avoid partial matches.
collapse_adjacent_fig_refs <- function(txt) {
  .run_track <- function(txt, singular_re, plural_str, singular_str) {
    ref_pat <- sprintf("\\*\\*%s \\d+\\*\\*", singular_re)
    run_pat <- sprintf("(%s)(?:\\s*;\\s*%s)+", ref_pat, ref_pat)
    repeat {
      m <- regexpr(run_pat, txt, perl = TRUE)
      if (m[1L] == -1L) break
      span  <- regmatches(txt, m)[[1L]]
      nums  <- sort(as.integer(regmatches(
                  span, gregexpr("[0-9]+", span, perl = TRUE))[[1L]]))
      label <- if (length(nums) > 1L) plural_str else singular_str
      repl  <- sprintf("**%s %s**", label, collapse_ints(nums))
      start <- m[1L]; end <- start + attr(m, "match.length") - 1L
      txt   <- paste0(substr(txt, 1L, start - 1L), repl,
                      substr(txt, end + 1L, nchar(txt)))
    }
    txt
  }
  txt <- .run_track(txt, "Extended Data Fig\\.",
                    "Extended Data Figs.", "Extended Data Fig.")
  txt <- .run_track(txt, "Fig\\.", "Figs.", "Fig.")
  txt
}

# ---- disk-mirror variant ---------------------------------------------------
# Walks the full sections/ directory, resolves narrative files in
# `section_order`, and writes a resolved mirror at sections/_resolved/.
# Files outside section_order (e.g. _setup.qmd, 06_figure_legends.qmd) are
# copied verbatim so Quarto includes from _resolved/ resolve correctly
# whether or not they contain tokens. Runs idempotently — safe to re-run.
resolve_sections_to_disk <- function(sections_dir    = NULL,
                                     out_dir         = NULL,
                                     section_order   = NULL,
                                     registry        = NULL,
                                     registry_path   = NULL) {
  if (is.null(sections_dir)) {
    sections_dir <- file.path(here::here(), "manuscript", "submission_prep", "sections")
  }
  if (is.null(out_dir)) {
    out_dir <- file.path(sections_dir, "_resolved")
  }
  if (is.null(section_order)) {
    # Default: the narrative order used by main_text.qmd.
    section_order <- c(
      "01_title_and_abstract",
      "02_introduction",
      "04_results",
      "05_discussion",
      "07_back_matter",
      "03_methods",
      "06_figure_legends"
    )
  }

  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Resolve the narrative files as a single pass so numbering is consistent
  # across the rendered document.
  narrative_paths <- file.path(sections_dir,
                               paste0(section_order, ".qmd"))
  narrative_paths <- narrative_paths[file.exists(narrative_paths)]
  resolution <- list(text = setNames(character(0L), character(0L)))
  if (length(narrative_paths) > 0L) {
    resolution <- resolve_sections(narrative_paths,
                                   registry      = registry,
                                   registry_path = registry_path)
  }

  # Write every .qmd in sections/ to out_dir — resolved if narrative,
  # verbatim otherwise. Skip the out_dir itself.
  all_qmd <- list.files(sections_dir, pattern = "\\.qmd$", full.names = TRUE)
  all_qmd <- all_qmd[!startsWith(all_qmd, out_dir)]
  for (p in all_qmd) {
    dest <- file.path(out_dir, basename(p))
    if (p %in% names(resolution$text)) {
      writeLines(resolution$text[[p]], dest)
    } else {
      file.copy(p, dest, overwrite = TRUE)
    }
  }

  # Persist the assignment maps so producer scripts and caption functions
  # can look up figure numbers / panel letters without re-parsing the
  # narrative. Sourced lazily from panel_assignments.R to avoid a tight
  # coupling between this file and that file's CLI guard.
  if (length(narrative_paths) > 0L && !is.null(resolution$maps)) {
    assignments <- resolution$maps
    if (is.null(assignments$meta)) {
      assignments$meta <- list(
        computed_at = Sys.time(),
        n_sections  = length(narrative_paths),
        n_cites     = NA_integer_,
        n_tokens    = length(resolution$token_index %||% list())
      )
    }
    cache_path <- file.path(here::here(), "results", "stats",
                            "panel_assignments.rds")
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(assignments, cache_path)
  }

  invisible(list(
    out_dir         = out_dir,
    resolved_paths  = file.path(out_dir, basename(narrative_paths)),
    maps            = if (length(narrative_paths) > 0L) resolution$maps else NULL
  ))
}

# ---- Registry-vs-disk integrity check -------------------------------------
# Ensures the registry is the hard source of truth for figure numbering:
#   * Every @fig/@edfig token used in prose must be defined in the registry
#     (enforced implicitly by the resolver — unknown tokens throw).
#   * Every registered EDFig group should have an on-disk output directory
#     (candidate names: results/ExtendedData/ed_<token>/ OR
#     results/ExtendedData/ExtendedData<NN>_<suffix>/). Missing dirs emit a
#     warning (soft gate) so the render keeps going during the disk refactor.
#   * Orphan disk dirs without any registry binding also emit a warning.
#   * Collisions (two dirs resolving to the same token) emit a hard error
#     regardless of strictness because they break the resolver's numbering.
#
# `strict = TRUE` (or `BAVM_REGISTRY_STRICT=on`) promotes the soft warnings
# to errors, gating the render on a clean registry↔disk alignment. Default
# off for the transition period.
check_registry_disk_alignment <- function(registry     = NULL,
                                          ed_dir       = NULL,
                                          strict       = NULL) {
  if (is.null(registry)) registry <- load_panel_registry()
  if (is.null(ed_dir)) {
    ed_dir <- file.path(here::here(), "results", "ExtendedData")
  }
  if (is.null(strict)) {
    # Default strict now that the disk is aligned to the registry (every
    # registered ED token has a canonical ed_<token>/ dir, no orphans).
    # Opt out with BAVM_REGISTRY_STRICT=off for debugging pre-rename states.
    strict <- !identical(tolower(Sys.getenv("BAVM_REGISTRY_STRICT")), "off")
  }

  ed_tokens <- names(registry$figures)[vapply(
    registry$figures, function(g) identical(g$track, "EDFig"), logical(1L)
  )]
  disk_dirs <- if (dir.exists(ed_dir)) {
    list.dirs(ed_dir, full.names = FALSE, recursive = FALSE)
  } else character(0L)
  # Dirs prefixed `_retired_` are legacy content intentionally parked on
  # disk (no longer registered, kept for reference). They're excluded from
  # the orphan check so they don't trip the gate. The `numbered/` sibling
  # is also excluded — it holds the manuscript-numbered symlinks
  # (edNN_<token>) maintained by sync_panel_prefixes.R, which point back
  # to the canonical ed_<token> dirs at the parent level.
  disk_dirs <- disk_dirs[!startsWith(disk_dirs, "_retired_") &
                         disk_dirs != "numbered"]

  # Match each registered token to 0+ disk dirs by trying the two canonical
  # patterns and falling back to a loose suffix match.
  match_one <- function(tok) {
    slug <- sub("^ed_", "", tok)
    candidates <- c(
      paste0("ed_", slug),                             # canonical
      paste0("ExtendedData_", slug),                   # legacy slug-only
      grep(paste0("_", slug, "$"), disk_dirs, value = TRUE)  # numbered+slug
    )
    intersect(unique(candidates), disk_dirs)
  }
  token_to_dirs <- setNames(lapply(ed_tokens, match_one), ed_tokens)
  missing_tokens <- ed_tokens[lengths(token_to_dirs) == 0L]
  collided      <- ed_tokens[lengths(token_to_dirs) > 1L]
  claimed_dirs  <- unlist(token_to_dirs, use.names = FALSE)
  orphan_dirs   <- setdiff(disk_dirs, claimed_dirs)

  diag <- list(
    ed_tokens        = ed_tokens,
    token_to_dirs    = token_to_dirs,
    missing_tokens   = missing_tokens,
    collided_tokens  = collided,
    orphan_dirs      = orphan_dirs,
    strict           = strict
  )

  fmt_list <- function(x) if (length(x) == 0L) "—" else paste(x, collapse = ", ")

  if (length(collided) > 0L) {
    stop(sprintf(
      "[registry] %d token(s) resolve to multiple disk dirs: %s. Pick one canonical path and remove/rename the others.",
      length(collided), fmt_list(collided)
    ), call. = FALSE)
  }

  warn_or_stop <- if (strict) stop else warning
  if (length(missing_tokens) > 0L) {
    warn_or_stop(sprintf(
      "[registry] %d registered ED token(s) have no disk directory: %s",
      length(missing_tokens), fmt_list(missing_tokens)
    ), call. = FALSE)
  }
  if (length(orphan_dirs) > 0L) {
    warn_or_stop(sprintf(
      "[registry] %d disk dir(s) under ExtendedData/ match no registry token: %s",
      length(orphan_dirs), fmt_list(orphan_dirs)
    ), call. = FALSE)
  }

  message(sprintf(
    "[registry] %d ED tokens | %d dirs matched | %d orphan dir(s) | %d missing dir(s)%s",
    length(ed_tokens),
    sum(lengths(token_to_dirs) > 0L),
    length(orphan_dirs),
    length(missing_tokens),
    if (strict) " [STRICT]" else ""
  ))

  invisible(diag)
}

# ---- Citation-coverage gate ------------------------------------------------
# Every group registered with track = "Fig" / "EDFig" must be cited at least
# once via @fig[...] / @edfig[...] in the canonical narrative sections (the
# same set compute_panel_assignments() walks). An uncited group has no number
# assigned, so any caption fn that calls figure_number("<that_key>") falls back
# to NA and the rendered docx prints "Extended Data Fig. NA" — silently. This
# check fails fast at render setup so the author sees the broken token before
# the docx ships.
#
# `strict = TRUE` (default; opt out via BAVM_CITATION_STRICT=off) promotes the
# coverage failure to an error. `strict = FALSE` warns instead — useful while
# bulk-renaming registry keys mid-edit.
#
# Returns invisibly the diagnostic list (all keys, cited keys, uncited keys).
check_registry_citation_coverage <- function(assignments = NULL,
                                             registry    = NULL,
                                             strict      = NULL) {
  if (is.null(registry))    registry    <- load_panel_registry()
  if (is.null(assignments)) {
    # Lazy load — most callers (e.g. _setup.qmd) already have the cache; fall
    # back to recomputing if not present.
    p <- file.path(here::here(), "results", "stats", "panel_assignments.rds")
    assignments <- if (file.exists(p)) readRDS(p) else {
      # Avoid circular dep on panel_assignments.R (which sources this file).
      # Inline the minimal recompute here.
      ns <- new.env(parent = baseenv())
      sys.source(file.path(here::here(), "analysis", "pipeline",
                           "helpers",
                           "panel_assignments.R"), envir = ns)
      ns$compute_panel_assignments(registry = registry)
    }
  }
  if (is.null(strict)) {
    strict <- !identical(tolower(Sys.getenv("BAVM_CITATION_STRICT")), "off")
  }

  fig_groups   <- names(Filter(function(g) identical(g$track, "Fig"),
                               registry$figures))
  edfig_groups <- names(Filter(function(g) identical(g$track, "EDFig"),
                               registry$figures))

  cited_fig   <- names(assignments$group_number$Fig)   %||% character(0L)
  cited_edfig <- names(assignments$group_number$EDFig) %||% character(0L)

  uncited_fig   <- setdiff(fig_groups,   cited_fig)
  uncited_edfig <- setdiff(edfig_groups, cited_edfig)

  diag <- list(
    fig_groups       = fig_groups,
    edfig_groups     = edfig_groups,
    cited_fig        = cited_fig,
    cited_edfig      = cited_edfig,
    uncited_fig      = uncited_fig,
    uncited_edfig    = uncited_edfig,
    strict           = strict
  )

  if (length(uncited_fig) > 0L || length(uncited_edfig) > 0L) {
    bullets <- c(
      if (length(uncited_fig) > 0L)
        sprintf("  Fig:   %s",   paste(uncited_fig,   collapse = ", ")),
      if (length(uncited_edfig) > 0L)
        sprintf("  EDFig: %s", paste(uncited_edfig, collapse = ", "))
    )
    msg <- paste(
      "[citation-coverage] registered figure groups have no @-citation in narrative prose:",
      paste(bullets, collapse = "\n"),
      "Fix: cite a panel of each via @fig[<token>] / @edfig[<token>] in the relevant 0X_*.qmd,",
      "or remove the entry from analysis/pipeline/panel_registry.R if no longer needed.",
      sep = "\n"
    )
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  } else {
    message(sprintf(
      "[citation-coverage] %d Fig + %d EDFig groups all cited%s",
      length(fig_groups), length(edfig_groups),
      if (strict) " [STRICT]" else ""
    ))
  }

  invisible(diag)
}

# ---- Supplementary body-order gate -----------------------------------------
# Walks manuscript/submission_prep/supplementary.qmd, extracts every inline
# `r captions$<fn>(stats)` call in source order, maps each caption fn back
# to the primary registry key it documents (the first `figure_number("ed_*")`
# call inside the fn body in captions.R), and looks up the resolver-assigned
# ED number. Errors if the sequence is not monotonically ascending — which
# means an embed block is in the wrong slot (the rendered docx would print
# non-monotonic ED numbers like "ED 1, ED 8, ED 3...").
#
# Together with check_registry_citation_coverage() this gives the two
# guarantees the pipeline needs:
#   coverage  : every registered group is cited at least once in prose,
#               so figure_number() never returns NA into a caption.
#   body order: the rendered ED sequence in supplementary.qmd is monotonic,
#               so the TOC and body always agree on the manuscript order.
#
# `strict = TRUE` (default; opt out via BAVM_BODY_ORDER_STRICT=off) makes
# misordering a hard error. Returns invisibly the diagnostic data frame.
check_supplementary_body_order <- function(
    sup_path     = NULL,
    captions_path = NULL,
    assignments  = NULL,
    strict       = NULL) {

  if (is.null(sup_path)) {
    sup_path <- file.path(here::here(), "manuscript", "submission_prep", "supplementary.qmd")
  }
  if (is.null(captions_path)) {
    captions_path <- file.path(here::here(), "manuscript", "submission_prep", "captions.R")
  }
  if (is.null(assignments)) {
    p <- file.path(here::here(), "results", "stats", "panel_assignments.rds")
    assignments <- if (file.exists(p)) readRDS(p) else {
      ns <- new.env(parent = baseenv())
      sys.source(file.path(here::here(), "analysis", "pipeline",
                           "helpers",
                           "panel_assignments.R"), envir = ns)
      ns$compute_panel_assignments()
    }
  }
  if (is.null(strict)) {
    strict <- !identical(tolower(Sys.getenv("BAVM_BODY_ORDER_STRICT")), "off")
  }

  sup <- readLines(sup_path, warn = FALSE)
  # Skip lines inside fenced code chunks (the TOC chunk references
  # captions$... false-positively otherwise).
  in_chunk <- FALSE
  fn_seq <- character(0L)
  for (i in seq_along(sup)) {
    if (grepl("^```", sup[i])) { in_chunk <- !in_chunk; next }
    if (in_chunk) next
    m <- regmatches(sup[i], regexpr("captions\\$[A-Za-z0-9_]+", sup[i]))
    if (length(m) == 1L && nzchar(m)) {
      fn_seq <- c(fn_seq, sub("captions\\$", "", m))
    }
  }
  # Restrict to ED captions; main-Fig captions live in main_text.qmd, not here.
  # Extracted later by the per-fn primary-key lookup falling back to NA when
  # the fn does not call figure_number("ed_*").

  # Parse captions.R to map fn name -> first ed_ key it numbers.
  cap <- readLines(captions_path, warn = FALSE)
  fn_starts <- grep("^\\s*[A-Za-z_][A-Za-z0-9_]*\\s*=\\s*function", cap)
  fn_names  <- sub("^\\s*([A-Za-z_][A-Za-z0-9_]*)\\s*=\\s*function.*", "\\1",
                   cap[fn_starts])

  fn_key <- function(fn_name) {
    idx <- which(fn_names == fn_name)
    if (length(idx) == 0L) return(NA_character_)
    start <- fn_starts[idx[1L]]
    end   <- if (idx[1L] < length(fn_starts)) fn_starts[idx[1L] + 1L] - 1L
             else length(cap)
    body  <- paste(cap[start:end], collapse = "\n")
    m <- regmatches(body,
                    regexpr('figure_number\\("(ed_[A-Za-z0-9_]+)"\\)',
                            body, perl = TRUE))
    if (length(m) == 0L || !nzchar(m)) return(NA_character_)
    sub('.*figure_number\\("(ed_[A-Za-z0-9_]+)"\\).*', "\\1", m)
  }

  keys <- vapply(fn_seq, fn_key, character(1L))
  ed_num <- assignments$group_number$EDFig
  nums <- unname(ed_num[keys])

  diag <- data.frame(
    pos          = seq_along(fn_seq),
    caption_fn   = fn_seq,
    primary_key  = ifelse(is.na(keys), "", keys),
    ed_number    = nums,
    stringsAsFactors = FALSE
  )

  # Drop fn calls that don't number an ED group (e.g. main-Fig captions
  # accidentally pasted in supplementary.qmd; or fns that legitimately have
  # no ed_* primary, in which case they're spurious here and should error).
  ed_rows <- diag[!is.na(diag$ed_number), , drop = FALSE]
  if (nrow(ed_rows) == 0L) return(invisible(diag))

  diffs <- diff(ed_rows$ed_number)
  bad <- which(diffs <= 0L)

  if (length(bad) > 0L) {
    bullets <- vapply(bad, function(k) {
      sprintf("  pos %d (%s, ED %d) -> pos %d (%s, ED %d)",
              ed_rows$pos[k], ed_rows$caption_fn[k], ed_rows$ed_number[k],
              ed_rows$pos[k + 1L], ed_rows$caption_fn[k + 1L],
              ed_rows$ed_number[k + 1L])
    }, character(1L))
    msg <- paste(
      "[body-order] supplementary.qmd ED embed sequence is not monotonic ascending:",
      paste(bullets, collapse = "\n"),
      "Fix: move the lower-numbered embed block before the higher-numbered one,",
      "or check whether a caption fn is referencing the wrong figure_number() key.",
      sep = "\n"
    )
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  } else {
    message(sprintf(
      "[body-order] %d ED embeds in supplementary.qmd are monotonic ascending%s",
      nrow(ed_rows), if (strict) " [STRICT]" else ""
    ))
  }

  invisible(diag)
}

# ---- Citation-order check --------------------------------------------------
# Verifies three ascending-order invariants and emits one warning per
# violation. Does NOT abort the render (strict = FALSE by default;
# set BAVM_CITE_ORDER_STRICT=on in the environment to promote to errors).
#
# Invariant 1 — group-level non-interleaving:
#   Every first-citation of any token in Fig N must precede every
#   first-citation of any token in Fig N+1 (no interleaving across groups).
#   Same rule for Extended Data Fig and for main Table tracks.
#
# Invariant 2 — panel-letter order within each group:
#   The first citation of panel 'a' must precede panel 'b', etc.
#   Fill-in panels (assigned a letter by the registry fill-in pass but never
#   cited in prose) are excluded — they have no prose position.
#
# Invariant 3 — SuppTable ascending:
#   Supplementary Table N first-cited before Table N+1.
check_citation_order <- function(all_citations, maps, token_index,
                                 strict = NULL) {
  if (is.null(strict)) {
    strict <- identical(tolower(Sys.getenv("BAVM_CITE_ORDER_STRICT")), "on")
  }

  # Sort key: unique position within the flattened document.
  .key <- function(sec, pos) as.numeric(sec) * 1e12 + as.numeric(pos)

  # --- Build first-citation table -------------------------------------------
  seen <- character(0L)
  fc_type  <- character(0L)   # "group" | "panel" | "table"
  fc_track <- character(0L)
  fc_group <- character(0L)   # NA_character_ for table rows
  fc_token <- character(0L)   # NA_character_ for group rows
  fc_key   <- numeric(0L)

  for (ci in all_citations) {
    track <- .macro_to_track(ci$macro)

    if (track %in% c("Fig", "EDFig")) {
      info <- token_index[[ci$tokens[1L]]]
      grp  <- info$group

      gkey <- paste0("G:", track, ":", grp)
      if (!gkey %in% seen) {
        seen     <- c(seen, gkey)
        fc_type  <- c(fc_type,  "group")
        fc_track <- c(fc_track, track)
        fc_group <- c(fc_group, grp)
        fc_token <- c(fc_token, NA_character_)
        fc_key   <- c(fc_key,   .key(ci$section_ord, ci$pos))
      }

      if (info$kind == "panel") {
        for (tok in ci$tokens) {
          pkey <- paste0("P:", grp, ":", tok)
          if (!pkey %in% seen) {
            seen     <- c(seen, pkey)
            fc_type  <- c(fc_type,  "panel")
            fc_track <- c(fc_track, track)
            fc_group <- c(fc_group, grp)
            fc_token <- c(fc_token, tok)
            fc_key   <- c(fc_key,   .key(ci$section_ord, ci$pos))
          }
        }
      }
    } else {
      for (tok in ci$tokens) {
        tkey <- paste0("T:", track, ":", tok)
        if (!tkey %in% seen) {
          seen     <- c(seen, tkey)
          fc_type  <- c(fc_type,  "table")
          fc_track <- c(fc_track, track)
          fc_group <- c(fc_group, NA_character_)
          fc_token <- c(fc_token, tok)
          fc_key   <- c(fc_key,   .key(ci$section_ord, ci$pos))
        }
      }
    }
  }

  if (length(fc_type) == 0L) {
    message("[citation-order] no citations to check")
    return(invisible(character(0L)))
  }

  violations <- character(0L)

  # --- Invariant 1: group-level non-interleaving ----------------------------
  for (track in c("Fig", "EDFig")) {
    gn <- maps$group_number[[track]]
    if (length(gn) < 2L) next
    label <- if (track == "Fig") "Fig." else "Extended Data Fig."
    grps  <- names(sort(unlist(gn)))

    for (i in seq_len(length(grps) - 1L)) {
      gA <- grps[i];       nA <- gn[[gA]]
      gB <- grps[i + 1L];  nB <- gn[[gB]]

      # All first-cite keys belonging to group A (panels + group-name refs)
      idx_A <- which((fc_type == "panel" | fc_type == "group") &
                       fc_track == track & !is.na(fc_group) & fc_group == gA)
      idx_B <- which((fc_type == "panel" | fc_type == "group") &
                       fc_track == track & !is.na(fc_group) & fc_group == gB)

      if (length(idx_A) == 0L || length(idx_B) == 0L) next

      first_B <- min(fc_key[idx_B])
      late_A  <- idx_A[fc_key[idx_A] > first_B]

      if (length(late_A) > 0L) {
        toks <- ifelse(!is.na(fc_token[late_A]),
                       fc_token[late_A],
                       paste0("[", gA, "]"))
        violations <- c(violations, sprintf(
          "%s %d→%d: token(s) of %s %d first-cited after first token of %s %d: %s",
          label, nA, nB, label, nA, label, nB,
          paste(toks, collapse = ", ")
        ))
      }
    }
  }

  # --- Invariant 2: panel-letter order within each group --------------------
  for (grp in names(maps$panel_letter)) {
    lm <- maps$panel_letter[[grp]]
    if (length(lm) < 2L) next

    tok1  <- names(lm)[1L]
    track <- token_index[[tok1]]$track
    fig_n <- maps$group_number[[track]][[grp]]
    label <- if (track == "Fig") "Fig." else "Extended Data Fig."

    toks_by_letter <- names(sort(unlist(lm)))  # a, b, c, ... order

    idx_grp   <- which(fc_type == "panel" & !is.na(fc_group) & fc_group == grp)
    prose_tok <- fc_token[idx_grp]
    prose_key <- fc_key[idx_grp]

    for (i in seq_len(length(toks_by_letter) - 1L)) {
      ta <- toks_by_letter[i]; tb <- toks_by_letter[i + 1L]
      hit_a <- which(prose_tok == ta); hit_b <- which(prose_tok == tb)
      if (length(hit_a) == 0L || length(hit_b) == 0L) next  # fill-in panel

      ka <- prose_key[hit_a[1L]]; kb <- prose_key[hit_b[1L]]
      if (ka > kb) {
        violations <- c(violations, sprintf(
          "%s %d: panel %s (%s) first-cited after panel %s (%s)",
          label, fig_n, lm[[ta]], ta, lm[[tb]], tb
        ))
      }
    }
  }

  # --- Invariant 3: SuppTable ascending -------------------------------------
  sn <- maps$table_number$SuppTable
  if (length(sn) >= 2L) {
    toks_by_num <- names(sort(unlist(sn)))
    idx_st  <- which(fc_type == "table" & fc_track == "SuppTable")
    st_tok  <- fc_token[idx_st]
    st_key  <- fc_key[idx_st]

    for (i in seq_len(length(toks_by_num) - 1L)) {
      ta <- toks_by_num[i]; tb <- toks_by_num[i + 1L]
      hit_a <- which(st_tok == ta); hit_b <- which(st_tok == tb)
      if (length(hit_a) == 0L || length(hit_b) == 0L) next

      ka <- st_key[hit_a[1L]]; kb <- st_key[hit_b[1L]]
      if (ka > kb) {
        violations <- c(violations, sprintf(
          "Supplementary Table %d (%s) first-cited after Table %d (%s)",
          sn[[ta]], ta, sn[[tb]], tb
        ))
      }
    }
  }

  # --- Report ---------------------------------------------------------------
  if (length(violations) > 0L) {
    msg <- paste(
      c(sprintf("[citation-order] %d violation(s):", length(violations)),
        paste0("  ", violations)),
      collapse = "\n"
    )
    if (strict) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  } else {
    message(sprintf(
      "[citation-order] OK — %d groups/panels/tables all cited in ascending order",
      length(fc_type)
    ))
  }

  invisible(violations)
}

# ---- CLI entry point -------------------------------------------------------
# `Rscript analysis/pipeline/helpers/resolve_panels.R` regenerates the mirror standalone.
# The guard detects direct-script invocation (only --file=... matches our own
# filename) so source()-ing this file from run_all.R does not re-trigger it.
.this_script <- "resolve_panels.R"
.invoked_as_script <- any(grepl(
  paste0("--file=.*/", .this_script, "$|--file=", .this_script, "$"),
  commandArgs(trailingOnly = FALSE)
))
if (.invoked_as_script) {
  info <- resolve_sections_to_disk()
  message(sprintf("✓ resolved %d narrative section(s) -> %s",
                  length(info$resolved_paths), info$out_dir))
}

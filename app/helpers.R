# app/helpers.R
# ---------------------------
# Constants and helpers for the Shiny app: drug colours, formatting, parameter
# sections, the verdict text, the state diagram, and input checks.
# ---------------------------

# ── Drugs ──────────────────────────────────────────────────────────────────────
# One fixed colour per drug, in a colour-blind-checked order. Charts list the drugs
# in this order so neighbouring series always get well-separated colours.
drug_cols <- c(Olanzapine     = "#2a78d6",
               Risperidone    = "#eb6834",
               Haloperidol    = "#1baf7a",
               Chlorpromazine = "#eda100",
               Quetiapine     = "#e87ba4")
drug_order <- names(drug_cols)
drug_code  <- c(Chlorpromazine = "CPZ", Haloperidol = "HAL", Quetiapine = "QUE",
                Risperidone = "RIS", Olanzapine = "OLA")
# second-line drug after each first-line drug (as in RunModel)
second_line <- c(Chlorpromazine = "Risperidone", Haloperidol = "Risperidone",
                 Quetiapine = "Haloperidol", Risperidone = "Haloperidol",
                 Olanzapine = "Haloperidol")
third_line_col <- "#7d8a90"   # fluphenazine decanoate, the third-line drug for every patient

ui_ink <- list(primary = "#1b2a33", secondary = "#4e5b63", muted = "#6f7a7f",
               grid = "#e3e7e3", axis = "#c5ccc8", surface = "#fcfcfb")
app_font <- "Atkinson Hyperlegible Next, system-ui, -apple-system, Segoe UI, sans-serif"

drug_tag <- function(drug) {
  tags$span(class = "drug", tags$span(class = "drug-dot", style = paste0("background:", drug_cols[[drug]])), drug)
}

# ── Formatting ─────────────────────────────────────────────────────────────────
minus <- "−"
fmt_usd <- function(x, digits = 0) {
  vapply(x, function(v) {
    if (is.na(v) || !is.finite(v)) return("–")
    s <- formatC(abs(v), format = "f", digits = digits, big.mark = ",")
    paste0(if (v < 0 && round(abs(v), digits) > 0) minus else "", "$", s)
  }, character(1))
}
fmt_num <- function(x, digits = 3) {
  vapply(x, function(v) {
    if (is.na(v) || !is.finite(v)) return("–")
    s <- formatC(abs(v), format = "f", digits = digits, big.mark = ",")
    paste0(if (v < 0 && round(abs(v), digits) > 0) minus else "", s)
  }, character(1))
}
fmt_pct <- function(x, digits = 0) paste0(formatC(100 * x, format = "f", digits = digits), "%")
# compact display of an input value (4 significant figures)
fmt_val <- function(x) {
  vapply(x, function(v) {
    if (is.na(v)) return("")
    if (v != 0 && abs(v) < 1e-3) return(formatC(v, format = "g", digits = 3))
    format(signif(v, 4), big.mark = ",", scientific = FALSE, trim = TRUE, drop0trailing = TRUE)
  }, character(1))
}
plural <- function(n, one, many = paste0(one, "s")) paste(format(n, big.mark = ","), if (n == 1) one else many)

# ── Parameter sections ─────────────────────────────────────────────────────────
# Sections follow the header rows of the Parameters sheet; blocks without a
# header row are named here.
param_sections <- function(p) {
  block_names <- c(pHospUntrt = "Hospitalisation and remission", age = "Model settings",
                   acuteDALY = "Disability weights", pEPSPLA = "Extrapyramidal symptoms (EPS)")
  header_names <- c(
    "Recover/Remission" = "Hospitalisation and remission",
    "Weight gain (>7%) probabilities" = "Weight gain",
    "Diabetes" = "Diabetes",
    "Personnel costs inputs" = "Staff time, inpatient care",
    "Personnel inputs for inpatient services" = "Staff time, inpatient care",
    "Personnel inputs for outpatient services" = "Staff time, outpatient care",
    "Personnel wages" = "Salaries and working hours",
    "Drug unit costs (per dosage unit)" = "Drug prices (per tablet or ampoule)",
    "Dosage units consumed per day" = "Daily doses",
    "Unit laboratory costs" = "Laboratory tests",
    "test frequency in a 3 month cycle" = "Laboratory tests",
    "Capital and overhead costs" = "Capital and overheads",
    "Travel costs" = "Travel, productivity and caregivers",
    "GDP per capita" = "Travel, productivity and caregivers",
    "General parameters" = "Travel, productivity and caregivers",
    "Prevalence of diabetes complications" = "Diabetes complications",
    "Disability weights for complications" = "Diabetes complications")
  grp <- rep(NA_character_, nrow(p)); cur <- NA_character_
  for (i in seq_len(nrow(p))) {
    v <- p$variable[i]; lab <- p$Parameter[i]
    if (is.na(v)) {
      cur <- if (is.na(lab)) NA_character_ else if (!is.na(header_names[lab])) header_names[[lab]] else lab
    } else {
      if (is.na(cur)) cur <- if (!is.na(block_names[v])) block_names[[v]] else "Other inputs"
      grp[i] <- cur
    }
  }
  grp
}
section_order <- c("Model settings", "Hospitalisation and remission", "Disability weights",
                   "Extrapyramidal symptoms (EPS)", "Weight gain", "Diabetes", "Diabetes complications",
                   "Drug prices (per tablet or ampoule)", "Daily doses", "Laboratory tests",
                   "Staff time, inpatient care", "Staff time, outpatient care", "Salaries and working hours",
                   "Capital and overheads", "Travel, productivity and caregivers", "Other inputs")

# Inputs the model actually reads (everything else in the sheet is shown but has no effect)
model_vars <- unique(c(all.names(body(RunModel)), all.names(body(mtrace))))

dist_labels <- c(fixed = "Fixed", normal = "Normal", lognormal = "Log-normal", beta = "Beta")
dist_summary <- function(dist, base, low, high, sem, alpha, beta) {
  mapply(function(d, b, l, h, s, a, bb) {
    switch(d %||% "",
      fixed     = "Fixed",
      normal    = paste0("Normal, SE ", fmt_val(s)),
      lognormal = paste0("Log-normal, 95% range ", fmt_val(l), "–", fmt_val(h)),
      beta      = paste0("Beta(", fmt_val(round(a, 2)), ", ", fmt_val(round(bb, 2)), ")"),
      "–")
  }, dist, base, low, high, sem, alpha, beta, USE.NAMES = FALSE)
}
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

# 95% interval of the probabilistic draws for one parameter
draw_interval <- function(dist, base, low, high, sem, alpha, beta) {
  switch(dist,
    normal    = qnorm(c(0.025, 0.975), base, sem),
    lognormal = qlnorm(c(0.025, 0.975), log(base), (log(high) - log(low)) / 3.92),
    beta      = qbeta(c(0.025, 0.975), alpha, beta),
    c(base, base))
}
# beta parameters from a mean and standard error (method of moments)
beta_mom <- function(m, se) {
  k <- m * (1 - m) / se^2 - 1
  c(alpha = m * k, beta = (1 - m) * k)
}

# Editors show values to 6 significant figures; an untouched value keeps its exact stored value
shown_row <- function(row, cols = c("basecase", "low", "high", "sem")) {
  for (cl in cols) row[[cl]] <- signif(row[[cl]], 6)
  row
}
keep_orig <- function(new, orig) {
  if (is.null(new) || length(new) == 0 || is.na(new)) return(NA_real_)
  if (!is.na(orig) && isTRUE(abs(new - signif(orig, 6)) <= 1e-12 * max(1, abs(orig)))) orig else new
}

# ── Checks on an edited parameter ──────────────────────────────────────────────
check_param <- function(dist, base, low, high, se) {
  e <- character(0)
  num <- function(x) !is.null(x) && length(x) == 1 && !is.na(x) && is.finite(x)
  if (!num(base)) return("Enter a base-case value.")
  has_low <- num(low); has_high <- num(high)
  if (has_low != has_high) e <- c(e, "Enter both a low and a high value, or leave both empty.")
  if (dist != "fixed" && !(has_low && has_high)) e <- c(e, "Enter a low and a high value for the one-way analysis.")
  if (has_low && has_high && !(low <= base && base <= high)) e <- c(e, "The base case must lie between the low and high values.")
  if (dist %in% c("normal", "beta") && !(num(se) && se > 0)) e <- c(e, "Enter a standard error above zero.")
  if (dist == "lognormal") {
    if (base <= 0) e <- c(e, "A log-normal input must be above zero.")
    if (has_low && low <= 0) e <- c(e, "The low value of a log-normal input must be above zero.")
    if (has_low && has_high && low == high) e <- c(e, "The low and high values must differ for a log-normal input.")
  }
  if (dist == "beta") {
    if (!(base > 0 && base < 1)) e <- c(e, "A beta input must lie between 0 and 1.")
    if (has_low && has_high && (low < 0 || high > 1)) e <- c(e, "The range of a beta input must lie within 0 and 1.")
    if (base > 0 && base < 1 && num(se) && se > 0 && se^2 >= base * (1 - base))
      e <- c(e, sprintf("The standard error is too large for a mean of %s; it must be below %s.",
                        fmt_val(base), fmt_val(sqrt(base * (1 - base)))))
  }
  unique(e)
}

# ── Which inputs differ from the published ones ────────────────────────────────
same_num <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) <= 1e-12 * pmax(1, abs(b)))
same_chr <- function(a, b) (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
rows_edited <- function(cur, pub, num_cols, chr_cols = character(0)) {
  ok <- rep(TRUE, nrow(cur))
  for (cl in num_cols) ok <- ok & same_num(cur[[cl]], pub[[cl]])
  for (cl in chr_cols) ok <- ok & same_chr(cur[[cl]], pub[[cl]])
  !ok
}
input_edits <- function(cur, pub) {
  list(
    parameters  = rows_edited(cur$parameters, pub$parameters, c("basecase", "low", "high", "sem", "alpha", "beta"), "dist"),
    transitions = rows_edited(cur$transitions, pub$transitions, c("basecase", "low", "high", "alpha")),
    lifeTab     = rows_edited(cur$lifeTab, pub$lifeTab, "p.death"),
    schizoTab   = rows_edited(cur$schizoTab, pub$schizoTab, "p.death"),
    lifeExp     = rows_edited(cur$lifeExp, pub$lifeExp, "yrs"),
    cpi         = rows_edited(cur$cpi, pub$cpi, "Value"),
    horizon     = cur$nCycles != pub$nCycles
  )
}
n_edits <- function(ed) sum(vapply(ed, sum, numeric(1)))

# ── The verdict: the conclusion stated in words ────────────────────────────────
verdict_text <- function(model, wtp) {
  r  <- model$ce_result
  ia <- IncrementalAnalysis(model)
  nmb <- -r$outcomes * wtp - r$costs
  best <- r$intervention[which.max(nmb)]
  frontier <- ia$ce_lines$intervention
  ib <- match(best, r$intervention)
  low <- function(x) tolower(x)

  if (length(frontier) == 1) {
    others <- r[r$intervention != frontier, ]
    nxt <- others$intervention[which.min(others$costs)]
    dc <- others$costs[which.min(others$costs)] - r$costs[ib]
    dd <- others$outcomes[which.min(others$costs)] - r$outcomes[ib]
    title <- sprintf("%s costs the least and leaves the fewest DALYs, so it dominates the other four drugs.", best)
    detail <- sprintf("The next cheapest, %s, costs %s more per patient and leaves %s more DALYs. That gap is small enough for the probabilistic analysis to matter.",
                      low(nxt), fmt_usd(dc, 2), fmt_num(dd, 3))
  } else {
    inc <- ia$incrementals
    k <- match(best, inc$intervention)
    if (k == 1) {
      title <- sprintf("At %s per DALY averted, %s is the most cost-effective drug.", fmt_usd(wtp), low(best))
      detail <- sprintf("It costs the least. %s averts more DALYs, at %s per DALY averted, which is above the threshold.",
                        inc$intervention[2], fmt_usd(inc$iCER[2]))
    } else {
      title <- sprintf("At %s per DALY averted, %s is the most cost-effective drug.", fmt_usd(wtp), low(best))
      detail <- sprintf("Compared with %s, it costs %s more per patient and averts %s DALYs: %s per DALY averted.",
                        low(inc$intervention[k - 1]), fmt_usd(inc$delta_costs[k], 2),
                        fmt_num(inc$delta_outcomes[k], 3), fmt_usd(inc$iCER[k]))
    }
  }
  list(best = best, title = title, detail = detail)
}

# ── Years spent in each health state ───────────────────────────────────────────
state_names <- c("Residual phase, on first-line drug", "Residual phase, on second-line drug",
                 "Residual phase, on third-line drug", "Stopped first-line drug", "Stopped second-line drug",
                 "Relapse after stopping first-line drug", "Relapse after stopping second-line drug",
                 "Relapse on first-line drug", "Relapse on second-line drug")
years_in_states <- function(trace) colSums(trace)[1:9]   # undiscounted, per patient

# The state-transition diagram, with the years a patient spends in each state
state_diagram <- function(years, first) {
  second <- second_line[[first]]
  c1 <- drug_cols[[first]]; c2 <- drug_cols[[second]]; c3 <- third_line_col
  esc <- htmltools::htmlEscape
  yrs <- function(i) paste0(formatC(years[i], format = "f", digits = 1), " years")

  node <- function(x, y, w, h, label, value, col, extra = "") {
    sprintf(paste0('<g><rect x="%d" y="%d" width="%d" height="%d" rx="10" fill="#fcfcfb" stroke="#c5ccc8"/>',
                   '<rect x="%d" y="%d" width="5" height="%d" rx="2.5" fill="%s"/>',
                   '<text x="%d" y="%d" class="d-label">%s</text>',
                   '<text x="%d" y="%d" class="d-value">%s</text>%s</g>'),
            x, y, w, h, x + 8, y + 10, h - 20, col, x + 22, y + 24, esc(label), x + 22, y + 48, esc(value), extra)
  }
  arrow <- function(d) sprintf('<path d="%s" class="d-edge" marker-end="url(#d-arrow)"/>', d)
  lab <- function(x, y, text, anchor = "start") sprintf('<text x="%d" y="%d" class="d-edge-label" text-anchor="%s">%s</text>', x, y, anchor, esc(text))
  head <- function(x, text, col) sprintf('<circle cx="%d" cy="24" r="6" fill="%s"/><text x="%d" y="29" class="d-head">%s</text>', x + 6, col, x + 18, esc(text))

  col_block <- function(xs, col, ri, oi, ai, ni, labels) {
    paste0(
      node(xs, 60, 136, 60, "Stopped drug", yrs(oi), col),
      node(xs + 164, 60, 136, 60, "Relapse, off drug", yrs(ai), col),
      node(xs, 196, 300, 60, "Residual phase, on drug", yrs(ri), col,
           if (labels) sprintf('<text x="%d" y="%d" class="d-note" text-anchor="end">starting state</text>', xs + 288, 196 + 48) else ""),
      node(xs + 164, 332, 136, 60, "Relapse, on drug", yrs(ni), col),
      arrow(sprintf("M%d,196 V122", xs + 68)),                     # stops for other reasons
      arrow(sprintf("M%d,90 H%d", xs + 136, xs + 162)),           # relapse while off treatment
      arrow(sprintf("M%d,122 V194", xs + 232)),                    # remission, restarts drug
      arrow(sprintf("M%d,256 V330", xs + 232)),                    # relapse on treatment
      arrow(sprintf("M%d,362 H%d V258", xs + 300, xs + 400)),      # remission, switches drug
      arrow(sprintf("M%d,226 H%d", xs + 300, xs + 330)),          # switches drug
      lab(xs + 76, 162, "stops for other reasons"),
      lab(xs + 240, 162, "remission, restarts"),
      lab(xs + 240, 298, "relapse"),
      lab(xs + 304, 384, "remission, switches drug"),
      lab(xs + 316, 188, "switch", "middle")
    )
  }

  svg <- paste0(
    '<svg viewBox="0 0 1000 480" role="img" aria-labelledby="d-title" xmlns="http://www.w3.org/2000/svg">',
    sprintf('<title id="d-title">Health states for patients who start on %s: years spent in each state</title>', esc(tolower(first))),
    '<defs><marker id="d-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto-start-reverse">',
    '<path d="M0,0 L10,5 L0,10 z" fill="#6f7a7f"/></marker>',
    '<style>.d-label{font:400 13px "Atkinson Hyperlegible Next",system-ui,sans-serif;fill:#4e5b63}',
    '.d-value{font:700 19px "Atkinson Hyperlegible Next",system-ui,sans-serif;fill:#1b2a33}',
    '.d-head{font:600 14px "Atkinson Hyperlegible Next",system-ui,sans-serif;fill:#1b2a33}',
    '.d-note{font:italic 400 12px "Atkinson Hyperlegible Next",system-ui,sans-serif;fill:#6f7a7f}',
    '.d-edge{fill:none;stroke:#6f7a7f;stroke-width:1.5}',
    '.d-edge-label{font:400 12px "Atkinson Hyperlegible Next",system-ui,sans-serif;fill:#6f7a7f}</style></defs>',
    head(20, paste("First line:", first), c1),
    head(352, paste("Second line:", second), c2),
    head(684, "Third line: fluphenazine decanoate", c3),
    col_block(20, c1, 1, 4, 6, 8, TRUE),
    col_block(352, c2, 2, 5, 7, 9, FALSE),
    node(684, 196, 296, 60, "Residual phase, on drug", yrs(3), c3),
    '<rect x="20" y="420" width="960" height="44" rx="10" fill="#e9ede9"/>',
    '<text x="40" y="447" class="d-label">Death is possible from every state, and the risk rises with age.</text>',
    '</svg>')
  HTML(svg)
}

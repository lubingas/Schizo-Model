# app/editors.R
# ---------------------------
# Dialogs for editing inputs: one parameter, one drug's set of transition
# probabilities, or one value in the mortality, life expectancy or price tables.
# ---------------------------

# Parameter ────────────────────────────────────────────────────────────────────
param_modal <- function(row, section, edited, used) {
  modalDialog(
    title = div(div(row$Parameter),
                div(class = "modal-sub", tags$span(class = "var-code", row$variable), paste0("  in ", section))),
    if (!used) div(class = "note", "The model doesn't use this input, so changing it won't change any result."),
    div(class = "ed-section",
        selectInput("pe_dist", "Distribution in the probabilistic analysis", width = "100%",
                    choices = c("Fixed (the same in every draw)" = "fixed", "Normal" = "normal",
                                "Log-normal" = "lognormal", "Beta (for probabilities and proportions)" = "beta"),
                    selected = row$dist)),
    div(class = "ed-section",
        h3("Values"),
        div(class = "ed-grid",
            numericInput("pe_base", "Base case", value = row$basecase),
            numericInput("pe_low", "Low, for the one-way analysis", value = row$low),
            numericInput("pe_high", "High, for the one-way analysis", value = row$high)),
        conditionalPanel("input.pe_dist == 'normal' || input.pe_dist == 'beta'",
                         div(class = "ed-grid", style = "margin-top: 0.6rem;",
                             numericInput("pe_se", "Standard error", value = row$sem))),
        div(style = "margin-top: 0.6rem;",
            checkboxInput("pe_scale", "Move the low, high and standard error in proportion when the base case changes", TRUE))),
    uiOutput("pe_feedback"),
    footer = tagList(
      if (edited) actionButton("pe_reset", "Reset to the published value", class = "btn-outline-secondary"),
      modalButton("Cancel"),
      uiOutput("pe_save_ui", inline = TRUE)),
    size = "l", easyClose = TRUE)
}

# Transition probabilities for one drug ─────────────────────────────────────────
tp_labels <- c(pHosp = "Relapses and goes to hospital", disAE = "Switches drug because of side effects",
               disEff = "Switches drug because it doesn't work", disOther = "Stops treatment for other reasons",
               rem = "Stays stable on the drug")
tp_label <- function(variable) unname(tp_labels[sub("(CPZ|HAL|OLA|RIS|QUE)$", "", variable)])

tp_modal <- function(drug, rows, edited) {
  n <- nrow(rows)
  cells <- lapply(seq_len(n), function(k) {
    tagList(div(class = "tp-label", tp_label(rows$variable[k])),
            numericInput(paste0("tp_base_", k), NULL, rows$basecase[k], min = 0, max = 1, step = 0.01),
            numericInput(paste0("tp_low_", k), NULL, rows$low[k], min = 0, max = 1, step = 0.01),
            numericInput(paste0("tp_high_", k), NULL, rows$high[k], min = 0, max = 1, step = 0.01),
            numericInput(paste0("tp_alpha_", k), NULL, rows$alpha[k], min = 0, step = 1))
  })
  modalDialog(
    title = div(div(paste("Yearly transitions on", tolower(drug))),
                div(class = "modal-sub", "Probabilities for a person in the residual phase on this drug; each column must add up to 1")),
    div(class = "tp-grid",
        div(class = "tp-head", "Each year, the person"), div(class = "tp-head", "Base case"),
        div(class = "tp-head", "Low set"), div(class = "tp-head", "High set"),
        div(class = "tp-head", "Count (Dirichlet)"),
        cells,
        div(class = "tp-head", "Total"), uiOutput("tp_sum_base"), uiOutput("tp_sum_low"), uiOutput("tp_sum_high"), div()),
    p(class = "small-print", style = "margin-top: 0.8rem;",
      "The one-way analysis uses the low and high sets together. The probabilistic analysis draws the five probabilities ",
      "together from a Dirichlet distribution with these counts, the number of patients in each outcome in the source trials. ",
      actionLink("tp_from_counts", "Set the base case from the counts.")),
    uiOutput("tp_feedback"),
    footer = tagList(
      if (edited) actionButton("tp_reset", "Reset to the published values", class = "btn-outline-secondary"),
      modalButton("Cancel"),
      uiOutput("tp_save_ui", inline = TRUE)),
    size = "l", easyClose = TRUE)
}

check_tp <- function(v) {
  e <- character(0)
  if (anyNA(v$base) || anyNA(v$low) || anyNA(v$high) || anyNA(v$alpha)) return("Fill in every cell.")
  for (nm in c("base", "low", "high")) {
    lab <- c(base = "base-case", low = "low", high = "high")[[nm]]
    if (any(v[[nm]] < 0 | v[[nm]] > 1)) e <- c(e, sprintf("The %s probabilities must lie between 0 and 1.", lab))
    s <- sum(v[[nm]])
    if (abs(s - 1) > 0.001) e <- c(e, sprintf("The %s probabilities add up to %s; they must add up to 1.", lab, formatC(s, format = "f", digits = 4)))
  }
  if (any(v$alpha <= 0)) e <- c(e, "Every count must be above zero.")
  e
}

# One value in a small table ───────────────────────────────────────────────────
value_tables <- list(
  lifeTab   = list(col = "p.death", what = "Probability of death in a year, general population", lo = 0, hi = 1),
  schizoTab = list(col = "p.death", what = "Extra probability of death in a year with schizophrenia", lo = 0, hi = 1),
  lifeExp   = list(col = "yrs", what = "Remaining life expectancy, years", lo = 0, hi = Inf),
  cpi       = list(col = "Value", what = "Consumer price index", lo = 0, hi = Inf)
)

value_modal <- function(table, row, edited) {
  spec <- value_tables[[table]]
  sub <- if (table == "cpi") paste("Year", row$Year) else paste("From age", row$ages)
  modalDialog(
    title = div(div(spec$what), div(class = "modal-sub", sub)),
    numericInput("ve_value", "Value", value = row[[spec$col]], width = "220px"),
    uiOutput("ve_feedback"),
    footer = tagList(
      if (edited) actionButton("ve_reset", "Reset to the published value", class = "btn-outline-secondary"),
      modalButton("Cancel"),
      uiOutput("ve_save_ui", inline = TRUE)),
    size = "m", easyClose = TRUE)
}

check_value <- function(table, x) {
  spec <- value_tables[[table]]
  if (is.null(x) || is.na(x)) return("Enter a value.")
  if (x < spec$lo || x > spec$hi) {
    return(if (is.finite(spec$hi)) sprintf("Enter a value between %s and %s.", spec$lo, spec$hi) else "Enter a value of zero or more.")
  }
  if (table == "cpi" && x <= 0) return("Enter a value above zero.")
  character(0)
}

save_button <- function(id, label, ok) {
  b <- actionButton(id, label, class = "btn-primary")
  if (!ok) b <- tagAppendAttributes(b, disabled = NA)
  b
}

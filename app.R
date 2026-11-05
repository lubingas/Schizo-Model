# app.R
# ---------------------------
# Shiny app for the schizophrenia cost-effectiveness model.
# Shiny sources the files in R/ (packages, inputs, model) before this file.
# Run locally with shiny::runApp() from the project folder.
# ---------------------------

library(shiny)
library(bslib)
library(reactable)
library(writexl)

source("app/helpers.R", local = TRUE)
source("app/charts.R", local = TRUE)
source("app/editors.R", local = TRUE)
source("app/equations.R", local = TRUE)

published <- read_inputs("data/data_input.xlsx")
published_wtp <- 1641    # 3 x Uganda's 2012 GDP per capita (US$547)
sections <- param_sections(published$parameters)

run_base <- function(inputs) {
  RunModel(inputs,
           baseparms = getparms(inputs$parameters, value = basecase),
           basetransitions = getparms(inputs$transitions, value = basecase))
}
published_model <- run_base(published)

rt_theme <- reactableTheme(
  color = ui_ink$primary, backgroundColor = "transparent", borderColor = "#dde2de",
  highlightColor = "#eef1ee", cellPadding = "9px 10px",
  style = list(fontFamily = "'Atkinson Hyperlegible Next', system-ui, sans-serif", fontSize = "0.92rem"),
  headerStyle = list(fontWeight = 600, color = ui_ink$secondary, fontSize = "0.84rem", borderBottomColor = "#c5ccc8"),
  groupHeaderStyle = list(fontWeight = 600),
  searchInputStyle = list(width = "100%", maxWidth = "420px", borderColor = "#c5ccc8", borderRadius = "8px",
                          padding = "6px 10px", marginBottom = "10px", backgroundColor = "#fcfcfb"),
  rowSelectedStyle = list(backgroundColor = "#eef1ee"))

edit_button <- function(input_id, value, label) {
  tags$button(type = "button", class = "btn btn-sm btn-outline-secondary", style = "padding: 0.15rem 0.7rem;",
              `aria-label` = label,
              onclick = sprintf("event.stopPropagation(); Shiny.setInputValue('%s', %s, {priority: 'event'})", input_id, value),
              "Edit")
}
row_click_js <- function(input_id, key = "row_id") {
  JS(sprintf("function(rowInfo) {
    if (!rowInfo) return;
    if (rowInfo.aggregated) { rowInfo.toggleRowExpanded(); return; }
    Shiny.setInputValue('%s', rowInfo.original.%s, {priority: 'event'});
  }", input_id, key))
}

# ── UI ────────────────────────────────────────────────────────────────────────

theme <- bs_theme(version = 5, bg = "#f2f4f1", fg = "#1b2a33", primary = "#1b2a33", secondary = "#4e5b63",
                  base_font = font_collection("Atkinson Hyperlegible Next", "system-ui", "sans-serif"),
                  heading_font = font_collection("Atkinson Hyperlegible Next", "system-ui", "sans-serif"),
                  "font-size-base" = "0.9375rem")

overview_panel <- nav_panel(
  "Overview", value = "overview",
  div(class = "overview-intro",
      div(h1("Which antipsychotic should people with schizophrenia in Uganda start on?"),
          p(class = "lede",
            "This model follows people with schizophrenia from age 25 for the rest of their lives and compares five drugs ",
            "they could start on: chlorpromazine, haloperidol, quetiapine, risperidone and olanzapine. It adds up the costs ",
            "of care, in 2012 US dollars, and the healthy life lost, in disability-adjusted life years (DALYs). Fewer DALYs is better.")),
      div(p("Each year a person can stay stable on their drug, switch because of side effects or poor response, stop for ",
            "other reasons, relapse and go to hospital, or die. Switching moves them to a second drug and then to a monthly ",
            "fluphenazine injection."),
          div(class = "drug-picker",
              radioButtons("ov_drug", "Show the years spent in each state for people who start on", inline = TRUE,
                           choiceNames = lapply(drug_order, drug_tag), choiceValues = drug_order,
                           selected = "Olanzapine")))),
  div(class = "diagram-wrap",
      uiOutput("diagram"),
      uiOutput("diagram_caption")),
  layout_columns(
    col_widths = c(8, 4), class = "mt-4",
    card(card_header(h2("Years in each state, by first drug"),
                     p(class = "card-sub", "Average years per person, not discounted")),
         card_body(div(style = "overflow-x: auto;", uiOutput("years_table")))),
    div(h2("The study"),
        div(class = "citation",
            p("Lubinga SJ, Mutamba BB, Nganizi A, Babigumira JB. A cost-effectiveness analysis of antipsychotics for ",
              "treatment of schizophrenia in Uganda. ", tags$em("Applied Health Economics and Health Policy"),
              " 2015;13(5):493–506. ",
              tags$a(href = "https://doi.org/10.1007/s40258-015-0176-3", target = "_blank", rel = "noopener",
                     "doi:10.1007/s40258-015-0176-3"))),
        p(class = "small-print", style = "margin-top: 1rem;",
          "This app runs a corrected R version of the published Excel model, so its results differ a little from the paper. ",
          actionLink("goto_about", "See what changed."),
          " Every input can be edited on the Inputs page, and the results update straight away."))))

inputs_panel <- nav_panel(
  "Inputs", value = "inputs",
  div(class = "inputs-toolbar",
      div(h1("Model inputs", style = "margin-bottom: 0.2rem;"), uiOutput("edit_summary")),
      div(class = "actions",
          actionButton("reset_all", "Reset all inputs", class = "btn-outline-secondary"),
          downloadButton("dl_inputs", "Download inputs (.xlsx)", class = "btn-outline-secondary"),
          div(class = "upload-input",
              fileInput("ul_inputs", NULL, accept = ".xlsx", buttonLabel = "Upload inputs (.xlsx)")))),
  navset_underline(
    id = "inputs_tab",
    nav_panel("Parameters", value = "params",
              p(class = "tab-intro",
                "Select a row to edit its base case, its range for the one-way analysis, or its distribution for the ",
                "probabilistic analysis. Inputs marked as not used are in the original spreadsheet but don't affect the model."),
              reactableOutput("tbl_params")),
    nav_panel("Transition probabilities", value = "tps",
              p(class = "tab-intro",
                "Each year, a person in the residual phase on a first- or second-line drug either relapses, switches drug, ",
                "stops treatment or stays stable. Select a drug to edit its five probabilities and counts together. ",
                "The one-way range of each probability comes from its count."),
              reactableOutput("tbl_tps")),
    nav_panel("Mortality and life expectancy", value = "mort",
              p(class = "tab-intro",
                "The WHO life table for Uganda gives the risk of dying within each 5-year age group; the model converts it ",
                "to a yearly risk. The extra yearly risk of death with schizophrenia is added to it, and the life expectancy ",
                "counts the years of life lost at each death. Select a row to edit it."),
              layout_columns(col_widths = c(4, 4, 4),
                             card(card_header(h3("General population, risk of death within each 5-year age group")), reactableOutput("tbl_lifeTab")),
                             card(card_header(h3("Extra yearly risk with schizophrenia")), reactableOutput("tbl_schizoTab")),
                             card(card_header(h3("Remaining life expectancy (years)")), reactableOutput("tbl_lifeExp")))),
    nav_panel("Price index", value = "cpi",
              p(class = "tab-intro",
                "Uganda's consumer price index. The model uses it to bring 2001 laboratory prices to 2012; the other prices ",
                "are already in 2012 US dollars."),
              div(style = "max-width: 420px;", card(reactableOutput("tbl_cpi")))),
    nav_panel("Model settings", value = "settings",
              layout_columns(
                col_widths = c(6, 6), class = "mt-3",
                card(card_header(h3("Time horizon")),
                     card_body(sliderInput("set_horizon", "Years followed, from the age of onset", min = 5, max = 70,
                                           value = published$nCycles, step = 1, width = "100%"),
                               p(class = "small-print", "The published model follows people for 70 years, to age 95."))),
                card(card_header(h3("Threshold for cost-effectiveness")),
                     card_body(numericInput("set_wtp", "Willingness to pay per DALY averted (US$)", value = published_wtp,
                                            min = 0, step = 100, width = "260px"),
                               p(class = "small-print",
                                 "The default is three times Uganda's 2012 GDP per capita (US$547), the World Health Organization's ",
                                 "upper threshold at the time. It sets the headline, the net benefit column, the net benefit in the one-way analysis and the marked point on the curves."))))))
)

base_panel <- nav_panel(
  "Base case", value = "base",
  card(card_header(h2("Costs and DALYs per patient"), uiOutput("base_sub")),
       card_body(reactableOutput("tbl_base"))),
  layout_columns(
    col_widths = c(8, 4), class = "mt-2",
    card(card_header(h2("Cost-effectiveness plane"),
                     p(class = "card-sub", "Each drug compared with the one that leaves the most DALYs")),
         card_body(plotlyOutput("plane", height = "460px"))),
    div(h3("Reading the plane"),
        p("A drug further right averts more DALYs; a drug lower down costs less. The line joins the drugs on the efficiency ",
          "frontier: the options worth considering, because no other drug, or mix of two, gives more health for less money."),
        p("Filled points are on the frontier. Hollow points are dominated."),
        p(class = "small-print", "Hover over a point for its costs and DALYs. The menu at the top right of the chart saves it as an image.")))
)

oneway_panel <- nav_panel(
  "One-way analysis", value = "oneway",
  card(card_body(
    div(class = "controls",
        selectInput("ow_drug", "Drug", choices = drug_order, selected = "Olanzapine", width = "180px"),
        selectInput("ow_comp", "Compared with", choices = drug_order, selected = "Risperidone", width = "180px"),
        selectInput("ow_outcome", "Result", width = "220px",
                    choices = c("Incremental cost" = "cost", "DALYs averted" = "outcome", "Cost per DALY averted" = "ICER",
                                "Net monetary benefit" = "NMB")),
        div(class = "ctl-wide", sliderInput("ow_n", "Inputs shown", min = 5, max = 40, value = 15, step = 1, width = "220px")),
        actionButton("ow_run", "Run the one-way analysis", class = "btn-primary")),
    p(class = "small-print", style = "margin: 0.8rem 0 0;",
      "Each input moves to its low and then its high value while every other input stays at its base case. ",
      "Each transition probability moves between the 2.5% and 97.5% points of its distribution from the trial counts, ",
      "and that drug's other transition probabilities are rescaled so they still add up to 1. ",
      "Net monetary benefit values each DALY averted at the willingness-to-pay threshold in the model settings and subtracts the extra cost. ",
      "A cost per DALY averted only ranks runs where the drug costs more and averts more; runs where it becomes dominant, ",
      "dominated, or cheaper but less effective are marked at the edges of the chart instead. ",
      "After a run, change the drugs or the result to see them straight away."))),
  uiOutput("ow_stale"),
  card(card_header(uiOutput("ow_title")), card_body(uiOutput("ow_body")))
)

psa_panel <- nav_panel(
  "Probabilistic analysis", value = "psa",
  card(card_body(
    div(class = "controls",
        numericInput("psa_n", "Number of draws", value = 1000, min = 100, max = 10000, step = 100, width = "150px"),
        numericInput("psa_seed", "Random seed", value = 1234, step = 1, width = "130px"),
        selectInput("psa_ref", "Planes compared with", choices = drug_order, selected = "Chlorpromazine", width = "190px"),
        numericInput("psa_maxwtp", "Curves up to (US$ per DALY)", value = 5000, min = 500, step = 500, width = "200px"),
        actionButton("psa_run", "Run the probabilistic analysis", class = "btn-primary")),
    p(class = "small-print", style = "margin: 0.8rem 0 0;",
      "Every input with a distribution is drawn at random at the same time, and the model runs once per draw. ",
      "1,000 draws take about 20 seconds."))),
  uiOutput("psa_stale"),
  uiOutput("psa_body")
)

about_panel <- nav_panel(
  "About the model", value = "about",
  div(class = "about-grid",
      div(h1("About the model"),
          p("A Markov cohort model with ten health states and yearly cycles. People start in the residual phase of schizophrenia ",
            "on one of five antipsychotics. Each year they can relapse and be hospitalised, switch to a second drug because of side ",
            "effects or poor response, stop treatment, or die. After the second drug, they move to monthly injections of fluphenazine ",
            "decanoate. Side effects (extrapyramidal symptoms, weight gain and diabetes) add to costs and to disability."),
          p("Costs cover staff time, drugs, laboratory tests, hospital capital and overheads, travel, and the time patients and ",
            "caregivers lose from work, in 2012 US dollars. Outcomes are disability-adjusted life years. Both are discounted ",
            "at 3% a year. The paper describes every input and its source.")),
      div(h2("Cite and reuse"),
          div(class = "citation",
              p("Lubinga SJ, Mutamba BB, Nganizi A, Babigumira JB. A cost-effectiveness analysis of antipsychotics for treatment ",
                "of schizophrenia in Uganda. ", tags$em("Applied Health Economics and Health Policy"), " 2015;13(5):493–506. ",
                tags$a(href = "https://doi.org/10.1007/s40258-015-0176-3", target = "_blank", rel = "noopener",
                       "doi:10.1007/s40258-015-0176-3"))),
          p(style = "margin-top: 1rem;", "The model code, its inputs and the original Excel model are on ",
            tags$a(href = "https://github.com/lubingas/Schizo-Model", target = "_blank", rel = "noopener", "GitHub", .noWS = "outside"), "."))),
  div(class = "mt-4",
      h2("Model structure"),
      div(class = "diagram-wrap",
          HTML(paste(readLines("www/model_structure.svg", warn = FALSE), collapse = "\n")),   # same file as in the README
          p(class = "diagram-caption",
            "Each year people stay in their state (the loops) or move along one arrow. Everyone starts in the residual phase on their ",
            "first-line drug. The Overview tab shows how many years people spend in each state."))),
  card(class = "mt-4",
       card_header(h2("How costs are calculated"),
                   p(class = "card-sub", "The cost of a year in each state; the model's cost functions follow these equations term by term")),
       card_body(class = "equations", cost_equations())),
  card(class = "mt-4",
       card_header(h2("How this version differs from the published model"),
                   p(class = "card-sub", "The published results came from the Excel model. This app runs an R version with these corrections.")),
       card_body(tags$table(class = "diff-table",
         tags$thead(tags$tr(tags$th("Change"), tags$th("Excel model (2015)"), tags$th("This version"),
                            tags$th(class = "num", "Effect on cost per patient"))),
         tags$tbody(
           tags$tr(tags$td(tags$strong("Cost of hospital stays")),
                   tags$td("The costs of the four relapse states were multiplied by disability-weighted occupancy instead of occupancy, which cut hospital costs by about a quarter"),
                   tags$td("Costs are multiplied by the share of people in each state"), tags$td(class = "num", "+$310 to +$415")),
           tags$tr(tags$td(tags$strong("HIV test")),
                   tags$td("The test's frequency was attached to the syphilis (VDRL) test, so the HIV test was never costed"),
                   tags$td("Costed at each hospital stay"), tags$td(class = "num", "+$50 to +$68")),
           tags$tr(tags$td(tags$strong("Laboratory price index")),
                   tags$td("A mistyped 2012 price index (196.64 instead of 196.46)"),
                   tags$td("All price indices are read from the price index table"), tags$td(class = "num", "−$0.65")),
           tags$tr(tags$td(tags$strong("Drug price years")),
                   tags$td("Chlorpromazine and haloperidol prices (MSH International Drug Price Indicator Guide 2013, already in US$) were inflated with an index of 145.14, and benzhexol was inflated twice"),
                   tags$td("MSH prices are used as published; benzhexol, from the National Medical Stores' October 2010 catalogue, is inflated once"),
                   tags$td(class = "num", "−$6 to −$17")),
           tags$tr(tags$td(tags$strong("Haloperidol dose and quetiapine price")),
                   tags$td("One 10 mg haloperidol tablet a day at $0.0206; quetiapine at $0.388 per tablet"),
                   tags$td("Two 5 mg tablets at $0.022 and quetiapine at $0.188, as in the paper's Table 2"),
                   tags$td(class = "num", "+$9 to +$12; quetiapine −$276")),
           tags$tr(tags$td(tags$strong("Benzhexol for outpatients")),
                   tags$td("Outpatients with extrapyramidal symptoms were costed for a 60-day inpatient supply"),
                   tags$td("Costed for the duration of the symptoms (0.25 years)"), tags$td(class = "num", "+$0.4 to +$0.5")),
           tags$tr(tags$td(tags$strong("Staff time")),
                   tags$td("Ward rounds are counted per week, but only the psychiatrist's round was divided by 7 for a daily cost, and the daily drug round was added in hours rather than dollars. Outpatient consultations costed the clinical officer and nurse for a full hour"),
                   tags$td("Every weekly round is divided by 7, the drug round is costed at the nurse's pay, and the consultation length applies to everyone at it. Staff time per hospital stay falls from $123.25 to $30.63"),
                   tags$td(class = "num", "−$257 to −$292")),
           tags$tr(tags$td(tags$strong("Probability of cost-effectiveness")),
                   tags$td("Net benefit counted each DALY as a gain, so the acceptability curves favoured drugs with more DALYs"),
                   tags$td("DALYs count as a loss: net benefit is minus the threshold times DALYs, minus cost"),
                   tags$td(class = "num", "None; changes the acceptability curves")),
           tags$tr(tags$td(tags$strong("Background mortality")),
                   tags$td("The WHO life table's probabilities of dying within each 5-year age group were used as yearly probabilities, so people died about five times too fast: life expectancy from age 25 was 16.7 years, against 39 in the same table"),
                   tags$td("Converted to yearly probabilities. Life expectancy from age 25 is 33.3 years, including the extra mortality from schizophrenia"),
                   tags$td(class = "num", "+$2,619 to +$2,657; 5.5 fewer DALYs"))))),
       card_footer("Inpatient drug costs also follow the length of stay, as in the Excel model; this only matters in the sensitivity analyses.")),
  card(class = "mt-3",
       card_header(h2("Published results and this version"),
                   p(class = "card-sub", "Per patient; costs in 2012 US$")),
       card_body(uiOutput("compare_table")))
)

ui <- page_navbar(
  id = "nav",
  title = tags$span(class = "brand",
                    tags$span(class = "brand-title", "Antipsychotics in Uganda"),
                    tags$span(class = "brand-sub", "Schizophrenia cost-effectiveness model")),
  window_title = "Antipsychotics for schizophrenia in Uganda: a cost-effectiveness model",
  theme = theme, fillable = FALSE, lang = "en",
  header = tagList(
    tags$head(tags$link(rel = "stylesheet", href = "app.css"),
              tags$link(rel = "icon", type = "image/svg+xml", href = "favicon.svg"),
              tags$meta(name = "description", content = "Interactive cost-effectiveness model of five antipsychotics for schizophrenia in Uganda (Lubinga et al. 2015).")),
    uiOutput("verdict")),
  overview_panel, inputs_panel, base_panel, oneway_panel, psa_panel, about_panel,
  footer = div(class = "app-footer",
               p("Lubinga SJ, Mutamba BB, Nganizi A, Babigumira JB. Applied Health Economics and Health Policy 2015;13(5):493–506. ",
                 "Model code on ", tags$a(href = "https://github.com/lubingas/Schizo-Model", target = "_blank", rel = "noopener", "GitHub", .noWS = "outside"), "."))
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  rv <- reactiveValues(inputs = published, edit_row = NULL, edit_orig = NULL,
                       tp_rows = NULL, ve = NULL)

  wtp <- reactive({
    w <- input$set_wtp
    if (is.null(w) || is.na(w) || w < 0) published_wtp else w
  })
  edits <- reactive(input_edits(rv$inputs, published))
  n_ed  <- reactive(n_edits(edits()))

  model <- reactive(tryCatch(run_base(rv$inputs), error = function(e) e))
  model_ok <- reactive({
    m <- model()
    if (inherits(m, "error")) validate(need(FALSE, paste("The model can't run with these inputs:", conditionMessage(m))))
    m
  })

  # ── Verdict strip ──────────────────────────────────────────────────────────
  output$verdict <- renderUI({
    m <- model()
    status <- if (n_ed() == 0) tags$span("Using the published inputs") else
      tagList(tags$span(class = "chip chip-edited", plural(n_ed(), "edited input")),
              actionLink("goto_inputs", "Review"))
    if (inherits(m, "error")) {
      return(div(class = "verdict is-error", role = "status",
                 div(class = "verdict-main",
                     div(p(class = "verdict-title", "The model can't run with these inputs."),
                         p(class = "verdict-detail", style = "color: inherit;", conditionMessage(m)))),
                 div(class = "verdict-status", status)))
    }
    v <- verdict_text(m, wtp())
    div(class = "verdict", role = "status", `aria-live` = "polite",
        div(class = "verdict-main",
            tags$span(class = "verdict-swatch", style = paste0("background:", drug_cols[[v$best]])),
            div(p(class = "verdict-title", v$title), p(class = "verdict-detail", v$detail))),
        div(class = "verdict-status", status))
  })
  observeEvent(input$goto_inputs, nav_select("nav", "inputs"))
  observeEvent(input$goto_about, nav_select("nav", "about"))

  # ── Overview ───────────────────────────────────────────────────────────────
  output$diagram <- renderUI({
    m <- model_ok()
    state_diagram(years_in_states(m$traces[[input$ov_drug]]), input$ov_drug)
  })
  output$diagram_caption <- renderUI({
    p(class = "diagram-caption",
      sprintf("Average years a person spends in each state over the %s-year time horizon, not discounted. Years spent on the second- and third-line drugs count from the switch.",
              rv$inputs$nCycles))
  })
  output$years_table <- renderUI({
    m <- model_ok()
    yrs <- sapply(drug_order, function(d) years_in_states(m$traces[[d]]))
    rows <- list(
      "Residual phase, on a drug" = colSums(yrs[1:3, , drop = FALSE]),
      "Stopped treatment"         = colSums(yrs[4:5, , drop = FALSE]),
      "In hospital after a relapse" = colSums(yrs[6:9, , drop = FALSE]),
      "Life expectancy from onset" = colSums(yrs))
    sel <- input$ov_drug
    th <- lapply(drug_order, function(d) tags$th(class = "num", style = if (d == sel) "color: var(--ink);" else NULL, drug_tag(d)))
    tags$table(class = "diff-table years-table",
               tags$thead(tags$tr(tags$th(""), th)),
               tags$tbody(lapply(names(rows), function(nm) {
                 tags$tr(tags$td(if (nm == "Life expectancy from onset") tags$strong(nm) else nm),
                         lapply(drug_order, function(d) tags$td(class = "num",
                           style = if (d == sel) "font-weight: 700;" else NULL,
                           formatC(rows[[nm]][[d]], format = "f", digits = 1))))
               })))
  })

  # ── Inputs: summary, reset, download, upload ───────────────────────────────
  output$edit_summary <- renderUI({
    n <- n_ed()
    if (n == 0) span(class = "edit-summary", "Using the published inputs.") else
      span(class = "edit-summary is-edited", paste0(plural(n, "input"), " changed from the published values."))
  })

  set_controls <- function(x, w = published_wtp) {
    updateSliderInput(session, "set_horizon", value = x$nCycles)
    updateNumericInput(session, "set_wtp", value = w)
  }
  observeEvent(input$reset_all, {
    rv$inputs <- published
    set_controls(published)
    showNotification("All inputs reset to the published values.", type = "message")
  })
  observeEvent(input$set_horizon, {
    if (!identical(as.numeric(input$set_horizon), as.numeric(rv$inputs$nCycles))) {
      x <- rv$inputs; x$nCycles <- as.integer(input$set_horizon); rv$inputs <- x
    }
  }, ignoreInit = TRUE)

  output$dl_inputs <- downloadHandler(
    filename = function() paste0("schizo-model-inputs-", Sys.Date(), ".xlsx"),
    content = function(file) {
      x <- rv$inputs
      writexl::write_xlsx(list(
        Parameters = x$parameters, `Transition Probabilities` = x$transitions, CPI = x$cpi,
        `Uganda Life Tables` = x$lifeTab, `Schizo Death Probability` = x$schizoTab, `Life Expectancy` = x$lifeExp,
        Settings = data.frame(setting = c("time_horizon_years", "wtp_per_daly"), value = c(x$nCycles, wtp()))), path = file)
    })

  observeEvent(input$ul_inputs, {
    f <- input$ul_inputs$datapath
    fail <- function(msg) showNotification(msg, type = "error", duration = 10)
    new <- tryCatch(read_inputs(f), error = function(e) NULL)
    if (is.null(new)) return(fail("That file can't be read. It needs the same sheets as the downloaded inputs file."))
    for (tb in c("parameters", "transitions", "cpi", "lifeTab", "schizoTab", "lifeExp")) {
      if (!all(names(published[[tb]]) %in% names(new[[tb]])))
        return(fail(sprintf("The %s sheet is missing columns. Start from the downloaded inputs file.", tb)))
      new[[tb]] <- new[[tb]][, names(published[[tb]])]
    }
    if (!setequal(na.omit(new$parameters$variable), na.omit(published$parameters$variable)) ||
        !setequal(na.omit(new$transitions$variable), na.omit(published$transitions$variable)))
      return(fail("The file's list of parameters doesn't match the model's. Edit values only, without adding or removing rows."))
    st <- tryCatch(readxl::read_excel(f, sheet = "Settings"), error = function(e) NULL)
    w <- wtp()
    if (!is.null(st)) {
      h <- st$value[st$setting == "time_horizon_years"]; if (length(h) == 1 && is.finite(h)) new$nCycles <- as.integer(min(max(h, 5), 70))
      ww <- st$value[st$setting == "wtp_per_daly"]; if (length(ww) == 1 && is.finite(ww)) w <- ww
    }
    test <- tryCatch(run_base(new), error = function(e) e)
    if (inherits(test, "error")) return(fail(paste("The model can't run with that file:", conditionMessage(test))))
    rv$inputs <- new
    set_controls(new, w)
    showNotification(sprintf("Loaded %s. %s differ from the published values.", input$ul_inputs$name,
                             plural(n_edits(input_edits(new, published)), "input")), type = "message")
  })

  # ── Inputs: parameters table and editor ────────────────────────────────────
  output$tbl_params <- renderReactable({
    p <- rv$inputs$parameters
    ed <- edits()$parameters
    keep <- !is.na(p$variable)
    d <- data.frame(row_id = which(keep), section = factor(sections[keep], levels = section_order),
                    name = p$Parameter[keep], code = p$variable[keep],
                    basecase = p$basecase[keep], low = p$low[keep], high = p$high[keep],
                    dist = dist_summary(p$dist[keep], p$basecase[keep], p$low[keep], p$high[keep], p$sem[keep], p$alpha[keep], p$beta[keep]),
                    used = p$variable[keep] %in% model_vars, edited = ed[keep], stringsAsFactors = FALSE)
    d <- d[order(d$section, d$row_id), ]
    reactable(
      d, theme = rt_theme, groupBy = "section", defaultExpanded = TRUE, searchable = TRUE, pagination = FALSE,
      highlight = TRUE, compact = FALSE, height = 680, onClick = row_click_js("param_click"),
      language = reactableLang(searchPlaceholder = "Search inputs by name or code"),
      rowClass = function(i) if (isTRUE(d$edited[i])) "row-edited" else NULL,
      columns = list(
        row_id = colDef(show = FALSE), code = colDef(show = FALSE), used = colDef(show = FALSE),
        section = colDef(name = "", width = 34, class = "group-col"),
        name = colDef(name = "Input", minWidth = 300, cell = function(value, i) {
          div(class = "param-name", div(value),
              div(class = "param-meta", tags$span(class = "var-code", d$code[i]),
                  if (!d$used[i]) tags$span(class = "chip chip-unused", "not used"),
                  if (isTRUE(d$edited[i])) tags$span(class = "chip chip-edited", "edited")))
        }),
        basecase = colDef(name = "Base case", align = "right", minWidth = 90, cell = function(v) fmt_val(v),
                          class = "num"),
        low = colDef(name = "Low", align = "right", minWidth = 80, cell = function(v) fmt_val(v), class = "num"),
        high = colDef(name = "High", align = "right", minWidth = 80, cell = function(v) fmt_val(v), class = "num"),
        dist = colDef(name = "Probabilistic analysis", minWidth = 190, style = list(color = ui_ink$secondary)),
        edited = colDef(name = "", width = 80, sortable = FALSE, align = "right",
                        cell = function(v, i) edit_button("param_click", d$row_id[i], paste("Edit", d$name[i]))))
    )
  })

  observeEvent(input$param_click, {
    i <- input$param_click
    row <- rv$inputs$parameters[i, ]
    rv$edit_row <- i; rv$edit_orig <- row
    showModal(param_modal(shown_row(row), sections[i], edited = edits()$parameters[i], used = row$variable %in% model_vars))
  })
  observeEvent(input$pe_base, {
    o <- rv$edit_orig; b <- input$pe_base
    if (is.null(o) || !isTRUE(input$pe_scale) || is.na(b) || is.na(o$basecase) || o$basecase == 0) return()
    k <- keep_orig(b, o$basecase) / o$basecase
    if (abs(k - 1) < 1e-12) return()
    if (!is.na(o$low))  updateNumericInput(session, "pe_low",  value = signif(o$low * k, 6))
    if (!is.na(o$high)) updateNumericInput(session, "pe_high", value = signif(o$high * k, 6))
    if (!is.na(o$sem))  updateNumericInput(session, "pe_se",   value = signif(o$sem * k, 6))
  })
  pe_vals <- reactive({
    o <- rv$edit_orig; req(o)
    list(dist = input$pe_dist, base = keep_orig(input$pe_base, o$basecase), low = keep_orig(input$pe_low, o$low),
         high = keep_orig(input$pe_high, o$high), se = keep_orig(input$pe_se, o$sem))
  })
  pe_errors <- reactive({
    v <- pe_vals(); req(v$dist)
    check_param(v$dist, v$base, v$low, v$high, v$se)
  })
  pe_ab <- reactive({     # beta parameters: keep the published counts unless the mean or SE changed
    v <- pe_vals(); o <- rv$edit_orig
    if (identical(o$dist, "beta") && isTRUE(all.equal(v$base, o$basecase)) && isTRUE(all.equal(v$se, o$sem)) && !is.na(o$alpha))
      c(alpha = o$alpha, beta = o$beta) else beta_mom(v$base, v$se)
  })
  output$pe_feedback <- renderUI({
    e <- pe_errors()
    if (length(e)) return(div(class = "ed-errors", tags$ul(lapply(e, tags$li))))
    v <- pe_vals()
    has_range <- !is.na(v$low %||% NA) && !is.na(v$high %||% NA)
    ab <- if (v$dist == "beta") pe_ab() else c(alpha = NA, beta = NA)
    iv <- draw_interval(v$dist, v$base, v$low, v$high, v$se, ab[["alpha"]], ab[["beta"]])
    txt <- switch(v$dist,
      fixed = "The probabilistic analysis uses the base case in every draw.",
      normal = sprintf("95%% of the probabilistic draws fall between %s and %s.", fmt_val(iv[1]), fmt_val(iv[2])),
      lognormal = sprintf("95%% of the probabilistic draws fall between %s and %s, the low and high values.", fmt_val(iv[1]), fmt_val(iv[2])),
      beta = sprintf("A beta distribution with α = %s and β = %s; 95%% of the probabilistic draws fall between %s and %s.",
                     fmt_val(round(ab[["alpha"]], 2)), fmt_val(round(ab[["beta"]], 2)), fmt_val(iv[1]), fmt_val(iv[2])))
    if (!has_range) txt <- paste(txt, "With no low and high values, the one-way analysis leaves this input out.")
    div(class = "ed-derived", txt)
  })
  output$pe_save_ui <- renderUI(save_button("pe_save", "Save changes", length(pe_errors()) == 0))
  observeEvent(input$pe_save, {
    if (length(pe_errors())) return()
    v <- pe_vals(); i <- rv$edit_row; o <- rv$edit_orig; x <- rv$inputs
    nz <- function(z) if (is.null(z) || is.na(z)) NA_real_ else z
    x$parameters$dist[i] <- v$dist
    x$parameters$basecase[i] <- v$base
    x$parameters$low[i] <- nz(v$low)
    x$parameters$high[i] <- nz(v$high)
    if (v$dist %in% c("normal", "beta")) x$parameters$sem[i] <- v$se
    if (v$dist == "beta") {
      ab <- pe_ab(); x$parameters$alpha[i] <- ab[["alpha"]]; x$parameters$beta[i] <- ab[["beta"]]
    }
    if (v$dist == "lognormal" && !(identical(o$dist, "lognormal") && isTRUE(all.equal(c(v$base, v$low, v$high), c(o$basecase, o$low, o$high))))) {
      x$parameters$alpha[i] <- log(v$base); x$parameters$beta[i] <- (log(v$high) - log(v$low)) / 3.92
    }
    rv$inputs <- x
    removeModal()
    showNotification(sprintf("Saved %s. The results have been updated.", o$variable), type = "message")
  })
  observeEvent(input$pe_reset, {
    i <- rv$edit_row; x <- rv$inputs
    cols <- c("basecase", "low", "high", "sem", "dist", "alpha", "beta")
    x$parameters[i, cols] <- published$parameters[i, cols]
    rv$inputs <- x
    removeModal()
    showNotification(sprintf("%s reset to the published value.", x$parameters$variable[i]), type = "message")
  })

  # ── Inputs: transition probabilities ───────────────────────────────────────
  output$tbl_tps <- renderReactable({
    t <- rv$inputs$transitions
    ed <- edits()$transitions
    keep <- which(!is.na(t$variable) & t$dist %in% "dirichlet")
    code <- stringr::str_extract(t$variable[keep], "(CPZ|HAL|OLA|RIS|QUE)$")
    drug <- names(drug_code)[match(code, drug_code)]
    # one-way range of each probability, from the counts (see getTPbounds)
    ow_low <- unlist(unname(getTPbounds(t, 0.025))); ow_high <- unlist(unname(getTPbounds(t, 0.975)))
    d <- data.frame(row_id = keep, drug = factor(drug, levels = drug_order), what = tp_label(t$variable[keep]),
                    code = t$variable[keep], basecase = t$basecase[keep], low = unname(ow_low[t$variable[keep]]),
                    high = unname(ow_high[t$variable[keep]]),
                    alpha = t$alpha[keep], edited = ed[keep], drug_name = drug, stringsAsFactors = FALSE)
    d <- d[order(d$drug, d$row_id), ]
    reactable(
      d, theme = rt_theme, groupBy = "drug", defaultExpanded = TRUE, pagination = FALSE, highlight = TRUE,
      onClick = row_click_js("tp_click", "drug_name"),
      rowClass = function(i) if (isTRUE(d$edited[i])) "row-edited" else NULL,
      columns = list(
        row_id = colDef(show = FALSE), code = colDef(show = FALSE), drug_name = colDef(show = FALSE),
        drug = colDef(name = "", width = 34, class = "group-col"),
        what = colDef(name = "Each year, the person", minWidth = 280, cell = function(v, i)
          div(class = "param-name", div(v), div(class = "param-meta", tags$span(class = "var-code", d$code[i]),
                                                  if (isTRUE(d$edited[i])) tags$span(class = "chip chip-edited", "edited")))),
        basecase = colDef(name = "Base case", align = "right", cell = function(v) fmt_val(v), class = "num"),
        low = colDef(name = "One-way low", align = "right", cell = function(v) fmt_val(v), class = "num"),
        high = colDef(name = "One-way high", align = "right", cell = function(v) fmt_val(v), class = "num"),
        alpha = colDef(name = "Count", align = "right", cell = function(v) fmt_val(v), class = "num"),
        edited = colDef(name = "", width = 80, sortable = FALSE, align = "right",
                        cell = function(v, i) edit_button("tp_click", paste0("'", d$drug_name[i], "'"),
                                                          paste("Edit transition probabilities for", d$drug_name[i]))))
    )
  })
  observeEvent(input$tp_click, {
    drug <- input$tp_click
    t <- rv$inputs$transitions
    idx <- which(!is.na(t$variable) & t$dist %in% "dirichlet" & grepl(paste0(drug_code[[drug]], "$"), t$variable))
    rv$tp_rows <- idx; rv$tp_drug <- drug
    showModal(tp_modal(drug, shown_row(t[idx, ], c("basecase", "alpha")), edited = any(edits()$transitions[idx])))
  })
  tp_vals <- reactive({
    n <- length(rv$tp_rows); req(n > 0)
    o <- rv$inputs$transitions[rv$tp_rows, ]
    g <- function(pre, col) vapply(seq_len(n), function(k) {
      z <- input[[paste0(pre, k)]]; if (is.null(z)) NA_real_ else keep_orig(as.numeric(z), o[[col]][k]) }, numeric(1))
    list(base = g("tp_base_", "basecase"), alpha = g("tp_alpha_", "alpha"))
  })
  tp_sum_ui <- function(nm) renderUI({
    s <- sum(tp_vals()[[nm]])
    bad <- is.na(s) || abs(s - 1) > 0.001
    div(class = paste("tp-sum", if (bad) "bad"), if (is.na(s)) "–" else formatC(s, format = "f", digits = 4))
  })
  output$tp_sum_base <- tp_sum_ui("base")
  tp_errors <- reactive(check_tp(tp_vals()))
  output$tp_feedback <- renderUI({
    e <- tp_errors()
    if (length(e)) div(class = "ed-errors", tags$ul(lapply(e, tags$li)))
  })
  output$tp_save_ui <- renderUI(save_button("tp_save", "Save changes", length(tp_errors()) == 0))
  observeEvent(input$tp_from_counts, {
    a <- tp_vals()$alpha
    if (anyNA(a) || sum(a) <= 0) return()
    for (k in seq_along(a)) updateNumericInput(session, paste0("tp_base_", k), value = signif(a[k] / sum(a), 6))
  })
  observeEvent(input$tp_save, {
    if (length(tp_errors())) return()
    v <- tp_vals(); idx <- rv$tp_rows; x <- rv$inputs
    x$transitions$basecase[idx] <- v$base; x$transitions$alpha[idx] <- v$alpha
    rv$inputs <- x
    removeModal()
    showNotification(sprintf("Saved the transition probabilities for %s.", tolower(rv$tp_drug)), type = "message")
  })
  observeEvent(input$tp_reset, {
    idx <- rv$tp_rows; x <- rv$inputs
    x$transitions[idx, c("basecase", "low", "high", "alpha")] <- published$transitions[idx, c("basecase", "low", "high", "alpha")]
    rv$inputs <- x
    removeModal()
    showNotification(sprintf("Transition probabilities for %s reset to the published values.", tolower(rv$tp_drug)), type = "message")
  })

  # ── Inputs: mortality, life expectancy and price index ─────────────────────
  value_table <- function(tb) renderReactable({
    spec <- value_tables[[tb]]
    x <- rv$inputs[[tb]]; ed <- edits()[[tb]]
    key <- if (tb == "cpi") "Year" else "ages"
    d <- data.frame(row_id = seq_len(nrow(x)), key = x[[key]], value = x[[spec$col]], edited = ed)
    reactable(d, theme = rt_theme, pagination = FALSE, highlight = TRUE, compact = TRUE,
              height = if (nrow(d) > 12) 460 else "auto", onClick = row_click_js(paste0("ve_click_", tb)),
              rowClass = function(i) if (isTRUE(d$edited[i])) "row-edited" else NULL,
              columns = list(
                row_id = colDef(show = FALSE),
                key = colDef(name = if (tb == "cpi") "Year" else "From age", align = "left", class = "num"),
                value = colDef(name = if (tb == "cpi") "Index" else if (tb == "lifeExp") "Years" else "Probability",
                               align = "right", class = "num", cell = function(v) fmt_val(v)),
                edited = colDef(name = "", width = 80, align = "right", sortable = FALSE,
                                cell = function(v, i) edit_button(paste0("ve_click_", tb), d$row_id[i], "Edit value"))))
  })
  for (tb in names(value_tables)) local({
    t <- tb
    output[[paste0("tbl_", t)]] <- value_table(t)
    observeEvent(input[[paste0("ve_click_", t)]], {
      i <- input[[paste0("ve_click_", t)]]
      rv$ve <- list(table = t, row = i)
      showModal(value_modal(t, shown_row(rv$inputs[[t]][i, ], value_tables[[t]]$col), edited = edits()[[t]][i]))
    })
  })
  ve_value <- reactive({ req(rv$ve); keep_orig(input$ve_value, rv$inputs[[rv$ve$table]][[value_tables[[rv$ve$table]]$col]][rv$ve$row]) })
  ve_errors <- reactive({ req(rv$ve); check_value(rv$ve$table, ve_value()) })
  output$ve_feedback <- renderUI({ e <- ve_errors(); if (length(e)) div(class = "ed-errors", e) })
  output$ve_save_ui <- renderUI(save_button("ve_save", "Save changes", length(ve_errors()) == 0))
  observeEvent(input$ve_save, {
    if (length(ve_errors())) return()
    t <- rv$ve$table; i <- rv$ve$row; x <- rv$inputs
    x[[t]][[value_tables[[t]]$col]][i] <- ve_value()
    rv$inputs <- x
    removeModal()
    showNotification("Saved. The results have been updated.", type = "message")
  })
  observeEvent(input$ve_reset, {
    t <- rv$ve$table; i <- rv$ve$row; x <- rv$inputs; cl <- value_tables[[t]]$col
    x[[t]][[cl]][i] <- published[[t]][[cl]][i]
    rv$inputs <- x
    removeModal()
    showNotification("Reset to the published value.", type = "message")
  })

  # ── Base case ──────────────────────────────────────────────────────────────
  output$base_sub <- renderUI({
    pr <- getparms(rv$inputs$parameters, value = basecase)
    p(class = "card-sub", sprintf("Over %s years from age %s, discounted at %s a year (costs) and %s a year (DALYs). Costs in 2012 US$.",
                                  rv$inputs$nCycles, fmt_val(pr[["age"]]), fmt_pct(pr[["discCosts"]], 1), fmt_pct(pr[["discOutcomes"]], 1)))
  })
  output$tbl_base <- renderReactable({
    m <- model_ok(); ia <- IncrementalAnalysis(m); w <- wtp()
    r <- m$ce_result
    inc <- ia$incrementals
    nmb <- -r$outcomes * w - r$costs
    best <- r$intervention[which.max(nmb)]; ib <- which.max(nmb)
    status <- ia$ce_points$status[match(r$intervention, ia$ce_points$intervention)]
    d <- data.frame(drug = r$intervention, cost = r$costs, dalys = r$outcomes, status = status,
                    dcost = r$costs - r$costs[ib], ddaly = r$outcomes - r$outcomes[ib],
                    icer = inc$iCER[match(r$intervention, inc$intervention)],
                    vs = NA_character_,
                    nmb = nmb - max(nmb), stringsAsFactors = FALSE)
    d$vs <- vapply(d$drug, function(x) { k <- match(x, inc$intervention); if (is.na(k) || k == 1) NA_character_ else inc$intervention[k - 1] }, character(1))
    d <- d[order(d$cost), ]
    show_icer <- nrow(inc) > 1
    lw <- tolower(best)
    reactable(
      d, theme = rt_theme, pagination = FALSE, sortable = TRUE,
      columns = list(
        drug = colDef(name = "Drug", minWidth = 160, cell = function(v) drug_tag(v)),
        cost = colDef(name = "Cost", align = "right", class = "num", cell = function(v) fmt_usd(v, 2)),
        dalys = colDef(name = "DALYs", align = "right", class = "num", cell = function(v) fmt_num(v, 3)),
        status = colDef(name = "Status", minWidth = 150, cell = function(v)
          tags$span(class = paste("chip", if (v == "On frontier") "chip-frontier" else "chip-dominated"),
                    switch(v, "On frontier" = "On the frontier", v))),
        dcost = colDef(name = paste("Cost compared with", lw), align = "right", class = "num", minWidth = 150,
                       cell = function(v, i) if (d$drug[i] == best) "\u2013" else paste0(if (v > 0) "+" else "", fmt_usd(v, 2))),
        ddaly = colDef(name = paste("DALYs compared with", lw), align = "right", class = "num", minWidth = 150,
                       cell = function(v, i) if (d$drug[i] == best) "\u2013" else paste0(if (v > 0) "+" else "", fmt_num(v, 3))),
        icer = colDef(name = "Cost per DALY averted", align = "right", class = "num", minWidth = 160, show = show_icer,
                      cell = function(v, i) if (is.na(v) || is.na(d$vs[i])) "\u2013" else
                        div(fmt_usd(v), div(class = "muted", style = "font-size: 0.8rem;", paste("vs", tolower(d$vs[i]))))),
        vs = colDef(show = FALSE),
        nmb = colDef(name = sprintf("Net benefit vs %s, at %s/DALY", lw, fmt_usd(w)), align = "right", class = "num", minWidth = 180,
                     cell = function(v, i) if (d$drug[i] == best) "\u2013" else fmt_usd(v, 2)))
    )
  })
  output$plane <- renderPlotly(plane_chart(model_ok(), session$clientData$output_plane_width %||% 1000))

  # ── One-way analysis ───────────────────────────────────────────────────────
  # the runs keep the costs and DALYs of every drug, so changing the drugs or the result
  # only re-picks from them (owsaPair); only a change of inputs needs a new run
  ow_state <- reactiveVal(NULL)
  run_owsa <- function() {
    m <- model()
    if (inherits(m, "error")) return(showNotification("Fix the inputs first: the model can't run with them.", type = "error"))
    inp <- rv$inputs
    lb <- getparms(inp$parameters, value = low);  hb <- getparms(inp$parameters, value = high)
    lt <- getTPbounds(inp$transitions, 0.025); ht <- getTPbounds(inp$transitions, 0.975)
    n <- length(lb) + length(unlist(lt))
    res <- withProgress(message = "Running the one-way analysis", value = 0, {
      withCallingHandlers(
        owsa(m, low_base = lb, low_transitions = lt, high_base = hb, high_transitions = ht),
        message = function(cnd) { incProgress(1 / n); invokeRestart("muffleMessage") })
    })
    tp_low <- unlist(unname(lt)); tp_high <- unlist(unname(ht))
    tp_drug <- names(drug_code)[match(stringr::str_extract(names(tp_low), "(CPZ|HAL|OLA|RIS|QUE)$"), drug_code)]
    info <- rbind(
      data.frame(code = names(lb), label = inp$parameters$Parameter[match(names(lb), inp$parameters$variable)],
                 low_in = fmt_val(lb), high_in = fmt_val(hb[names(lb)]), stringsAsFactors = FALSE),
      data.frame(code = names(tp_low), label = paste0(tp_drug, ": ", tolower(tp_label(names(tp_low)))),
                 low_in = fmt_val(tp_low), high_in = fmt_val(tp_high[names(tp_low)]), stringsAsFactors = FALSE))
    ow_state(list(res = res, info = info, inputs = inp))
  }
  observeEvent(input$ow_run, run_owsa())

  # the chosen drugs and result, picked from the stored runs
  ow_pair <- reactive({
    s <- ow_state(); req(s)
    validate(need(input$ow_drug != input$ow_comp, "Choose two different drugs."))
    owsaPair(s$res, input$ow_outcome, input$ow_drug, input$ow_comp, wtp = wtp())
  })

  output$ow_stale <- renderUI({
    s <- ow_state(); if (is.null(s)) return(NULL)
    if (!identical(s$inputs, rv$inputs)) div(class = "note note-stale", "The inputs have changed since this run. Run the one-way analysis again to update it.")
  })
  output$ow_title <- renderUI({
    s <- ow_state(); if (is.null(s) || input$ow_drug == input$ow_comp) return(h2("One-way sensitivity analysis"))
    ow <- ow_pair()
    what <- switch(input$ow_outcome, cost = "incremental cost", outcome = "DALYs averted", ICER = "cost per DALY averted",
                   NMB = sprintf("net monetary benefit (at %s per DALY averted)", fmt_usd(wtp())))
    sub <- if (input$ow_outcome == "ICER" && !is.na(ow$owsa_basestatus)) {
      "At the base case there is no cost per DALY averted to vary; see below."
    } else {
      sprintf("The %s inputs with the largest effect, out of %s. The vertical line is the base case: %s.",
              min(input$ow_n, nrow(ow$owsa_data)), nrow(ow$owsa_data), owsa_value_fmt(input$ow_outcome)(ow$owsa_basecase))
    }
    tagList(h2(sprintf("What moves the %s of %s compared with %s", what, tolower(input$ow_drug), tolower(input$ow_comp))),
            p(class = "card-sub", sub))
  })
  output$ow_body <- renderUI({
    s <- ow_state()
    if (is.null(s)) return(div(class = "empty-state", p("Run the one-way analysis to see which inputs move the result most.")))
    if (input$ow_drug == input$ow_comp) return(div(class = "empty-state", p("Choose two different drugs.")))
    ow <- ow_pair()
    if (input$ow_outcome == "ICER" && !is.na(ow$owsa_basestatus)) {
      drug <- tolower(input$ow_drug); comp <- tolower(input$ow_comp)
      msg <- switch(ow$owsa_basestatus,
                    dominant = sprintf("At the base case %s costs less than %s and averts at least as many DALYs, so it is dominant and there is no cost per DALY averted.", drug, comp),
                    dominated = sprintf("At the base case %s costs more than %s and averts no more DALYs, so it is dominated and there is no cost per DALY averted.", drug, comp),
                    sprintf("At the base case %s costs less than %s but averts fewer DALYs. Swap the two drugs to see the cost per DALY averted of %s compared with %s.", drug, comp, comp, drug))
      return(div(class = "empty-state", p(msg), p("Net monetary benefit, incremental cost or DALYs averted show what moves this comparison.")))
    }
    n <- min(input$ow_n, length(s$res$variable))
    tagList(plotlyOutput("tornado", height = paste0(105 + 28 * n, "px")),
            tags$details(style = "margin-top: 0.8rem;", tags$summary("Show as a table"), reactableOutput("ow_table")))
  })
  output$tornado <- renderPlotly({
    s <- ow_state(); req(s)
    tornado_chart(ow_pair(), s$info, input$ow_n, input$ow_outcome, input$ow_drug, input$ow_comp,
                  session$clientData$output_tornado_width %||% 1000)
  })
  output$ow_table <- renderReactable({
    s <- ow_state(); req(s)
    d <- tornado_data(ow_pair(), s$info, input$ow_n); d <- d[rev(seq_len(nrow(d))), ]
    f <- owsa_value_fmt(input$ow_outcome)
    res <- function(x, st) ifelse(is.na(st), f(x), owsa_status_text[st])   # a status where there is no ratio
    reactable(data.frame(Input = d$label, `Low input` = d$low_in, `Result at low` = res(d$low, d$low_status),
                         `High input` = d$high_in, `Result at high` = res(d$high, d$high_status), check.names = FALSE),
              theme = rt_theme, pagination = FALSE, compact = TRUE, sortable = FALSE,
              columns = list(Input = colDef(minWidth = 280)),
              defaultColDef = colDef(align = "right", class = "num", minWidth = 110))
  })

  # ── Probabilistic analysis ─────────────────────────────────────────────────
  psa_state <- reactiveVal(NULL)
  run_psa <- function() {
    m <- model()
    if (inherits(m, "error")) return(showNotification("Fix the inputs first: the model can't run with them.", type = "error"))
    n <- input$psa_n; seed <- input$psa_seed
    if (is.na(n) || n < 100 || n > 10000) return(showNotification("Choose between 100 and 10,000 draws.", type = "warning"))
    if (is.na(seed)) seed <- 1234
    set.seed(seed)
    res <- tryCatch(withProgressShiny(message = "Running the probabilistic analysis", detail = "Starting …", value = 0,
                                      RunPSA(m, nsims = round(n))),
                    error = function(e) e)
    if (inherits(res, "error")) return(showNotification(paste("The probabilistic analysis stopped:", conditionMessage(res)), type = "error", duration = 10))
    psa_state(list(res = res, inputs = rv$inputs, n = round(n), seed = seed))
  }
  observeEvent(input$psa_run, run_psa())

  output$psa_stale <- renderUI({
    s <- psa_state(); if (is.null(s)) return(NULL)
    if (!identical(s$inputs, rv$inputs)) div(class = "note note-stale", "The inputs have changed since this run. Run the probabilistic analysis again to update it.")
  })
  ceac_data <- reactive({
    s <- psa_state(); req(s)
    mx <- input$psa_maxwtp; if (is.null(mx) || is.na(mx) || mx <= 0) mx <- 5000
    RunCEAC(s$res, wtp = mx, by = mx / 100)
  })
  output$psa_body <- renderUI({
    s <- psa_state()
    if (is.null(s)) return(card(card_body(div(class = "empty-state", p("Run the probabilistic analysis to see how sure the result is.")))))
    tagList(
      layout_columns(
        col_widths = c(7, 5),
        card(card_header(h2("Chance each drug is the most cost-effective"),
                         p(class = "card-sub", sprintf("%s draws, seed %s", format(s$n, big.mark = ","), s$seed))),
             card_body(plotlyOutput("ceac", height = "430px"))),
        card(card_header(h2(sprintf("At %s per DALY averted", fmt_usd(wtp()))),
                         p(class = "card-sub", "Means and 95% intervals across the draws")),
             card_body(reactableOutput("psa_summary")))),
      card(class = "mt-2",
           card_header(h2(sprintf("Each drug compared with %s, draw by draw", tolower(input$psa_ref))),
                       p(class = "card-sub", "Points below the dotted threshold line are cost-effective at the threshold; the large point is the mean")),
           card_body(plotlyOutput("psa_planes", height = "560px"))))
  })
  output$ceac <- renderPlotly(ceac_chart(ceac_data(), wtp()))
  output$psa_planes <- renderPlotly({ s <- psa_state(); req(s); psa_plane_chart(s$res, input$psa_ref, wtp()) })
  output$psa_summary <- renderReactable({
    s <- psa_state(); req(s)
    pc <- genProbCE(s$res, wtp())
    q <- function(x) sprintf("%s to %s", x[1], x[2])
    d <- do.call(rbind, lapply(drug_order, function(dn) {
      cst <- s$res$costs[[dn]]; dl <- s$res$outcomes[[dn]]
      data.frame(drug = dn, cost = mean(cst), cost_ci = q(fmt_usd(quantile(cst, c(0.025, 0.975)))),
                 dalys = mean(dl), dalys_ci = q(fmt_num(quantile(dl, c(0.025, 0.975)), 2)),
                 p = pc[[dn]], stringsAsFactors = FALSE)
    }))
    reactable(d, theme = rt_theme, pagination = FALSE, sortable = TRUE, compact = TRUE,
              defaultSorted = list(p = "desc"),
              columns = list(
                drug = colDef(name = "Drug", minWidth = 125, cell = function(v) drug_tag(v)),
                cost = colDef(name = "Cost", align = "right", class = "num", minWidth = 112,
                              cell = function(v, i) div(fmt_usd(v), div(class = "muted", style = "font-size: 0.8rem;", d$cost_ci[i]))),
                cost_ci = colDef(show = FALSE), dalys_ci = colDef(show = FALSE),
                dalys = colDef(name = "DALYs", align = "right", class = "num", minWidth = 100,
                               cell = function(v, i) div(fmt_num(v, 2), div(class = "muted", style = "font-size: 0.8rem;", d$dalys_ci[i]))),
                p = colDef(name = "Chance best", align = "right", class = "num", minWidth = 80,
                           cell = function(v) fmt_pct(v))))
  })

  # ── About: published results and this version ──────────────────────────────
  output$compare_table <- renderUI({
    excel <- c(Chlorpromazine = 4617, Haloperidol = 4440, Quetiapine = 4943, Risperidone = 4424, Olanzapine = 4454)
    excel_d <- c(Chlorpromazine = 27.608, Haloperidol = 27.563, Quetiapine = 27.562, Risperidone = 27.557, Olanzapine = 27.552)
    pub <- setNames(published_model$ce_result$costs, published_model$ce_result$intervention)
    pubd <- setNames(published_model$ce_result$outcomes, published_model$ce_result$intervention)
    m <- model()
    live_ok <- !inherits(m, "error") && n_ed() > 0
    live <- if (live_ok) setNames(m$ce_result$costs, m$ce_result$intervention) else NULL
    lived <- if (live_ok) setNames(m$ce_result$outcomes, m$ce_result$intervention) else NULL
    tags$table(class = "diff-table",
      tags$thead(tags$tr(tags$th("Drug"),
                         tags$th(class = "num", "Cost, paper"), tags$th(class = "num", "Cost, this version"),
                         if (live_ok) tags$th(class = "num", "Cost, your inputs"),
                         tags$th(class = "num", "DALYs, paper"), tags$th(class = "num", "DALYs, this version"),
                         if (live_ok) tags$th(class = "num", "DALYs, your inputs"))),
      tags$tbody(lapply(drug_order, function(d) tags$tr(
        tags$td(drug_tag(d)), tags$td(class = "num", fmt_usd(excel[[d]])), tags$td(class = "num", fmt_usd(pub[[d]])),
        if (live_ok) tags$td(class = "num", fmt_usd(live[[d]])),
        tags$td(class = "num", fmt_num(excel_d[[d]], 3)), tags$td(class = "num", fmt_num(pubd[[d]], 3)),
        if (live_ok) tags$td(class = "num", fmt_num(lived[[d]], 3))))))
  })
}

shinyApp(ui, server)

# R/data_input.R
# ---------------------------
# Reads all model inputs into one list
# ---------------------------

# inputs no equation reads (none were used in the Excel model either): remission on haloperidol
# (patients in hospital are on chlorpromazine, so pRemitCPZ applies to every drug), death in hospital,
# ischaemic heart disease as a diabetes complication, cohort size, the share of patients seen by the
# pharmacy assistant (1), and the exchange rate (prices are already in US$)
unused_parameters <- c("pRemitHAL", "pDieHosp", "pIHD", "dIHD",
                       "cohort", "pharm_prop", "exchangerate_dollar")

read_inputs <- function(path = here::here("data", "data_input.xlsx")) {
  required_sheets <- c(
    parameters = "Parameters",
    transitions = "Transition Probabilities",
    cpi = "CPI",
    lifeTab = "Uganda Life Tables",
    schizoTab = "Schizo Death Probability",
    lifeExp = "Life Expectancy"
  )

  inputs <- purrr::map(required_sheets, ~ readxl::read_excel(path, sheet = .x))

  # parameters in the spreadsheet that the model doesn't use: left out, so they are not drawn in the
  # PSA, varied in the one-way analysis or listed in the app (they stay in data_input.xlsx)
  inputs$parameters <- dplyr::filter(inputs$parameters, !variable %in% unused_parameters)

  # model settings; these are not varied in sensitivity analyses
  c(inputs,
    list(nStates = 10,
         nCycles = 70,
         interventions = c("Chlorpromazine", "Haloperidol", "Quetiapine", "Risperidone", "Olanzapine")))
}

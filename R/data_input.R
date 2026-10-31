# R/data_input.R
# ---------------------------
# Reads all model inputs into one list
# ---------------------------

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

  # model settings; these are not varied in sensitivity analyses
  c(inputs,
    list(nStates = 10,
         nCycles = 70,
         interventions = c("Chlorpromazine", "Haloperidol", "Quetiapine", "Risperidone", "Olanzapine")))
}

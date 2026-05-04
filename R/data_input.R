
# some preliminary data to run model
# these are not varied in sensitivity analyses
nStates <- 10
nCycles <- 70
interventions <- c("Chlorpromazine", "Haloperidol", "Quetiapine", "Risperidone", "Olanzapine")

read_inputs <- function(path = here::here("data", "data_input.xlsx")) {
  required_sheets <- c(
    parameters = "Parameters",
    transitions = "Transition Probabilities",
    cpi = "CPI",
    lifeTab = "Uganda Life Tables",
    schizoTab = "Schizo Death Probability",
    lifeExp = "Life Expectancy"
  )
  
  purrr::imap(required_sheets, ~ readxl::read_excel(path, sheet = .x))

}


# inputs for PSA


# colors
col <- brewer.pal(6,"Dark2")



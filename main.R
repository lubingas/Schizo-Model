# main.R
# ---------------------------
# Primary driver script for model
# ---------------------------

# clear memory
rm(list = ls())

# Reproducibility
set.seed(1234)

# Load packages and helper functions
source(here::here("R", "00_packages.R"))
source(here::here("R", "data_input.R"))
source(here::here("R", "functions.R"))

# load data and place in global environment for other functions to use
list2env(read_inputs(), envir = .GlobalEnv)

## run model
base_case <- RunModel(intervention = interventions, 
                      baseparms = getparms(parameters, value = basecase), 
                      basetransitions = getparms(transitions, value = basecase), 
                      DALYs = TRUE)

# perform incremental analysis (establish efficiency frontier)
CEAnalysis <- IncrementalAnalysis(model = base_case)

# plot incremental analysis
ScatterPlot(CEAnalysis)

## ONE WAY SENSITIVITY ANALYSES
# iterate over variables and perform OWSA
dsa_res <- owsa(model = base_case,
                low_base = getparms(parameters, value = low),
                low_transitions = getTPset(transitions, value = low),
                high_base = getparms(parameters, value = high),
                high_transitions = getTPset(transitions, value = high),
                outcome = "cost", comparator = "Haloperidol")

# plot OWSA
TornadoPlot(owsa_list = dsa_res, max_vars = 30)


## PSA
# how many simulations?
nsims <- 1000

# run PSA and CEAC
#profvis({
psa_res <- RunPSA(
  
  model          = base_case,
  nsims          = nsims,         # tibble of PSA draws
  interventions  = interventions,
  parameters     = parameters,
  transitions    = transitions
  )
#})

psa_res$costs    # nsims × n_interventions
psa_res$outcomes # nsims × n_interventions

plot_ly() %>%
  add_trace(x = ~psa_res$outcomes$Chlorpromazine, y = ~psa_res$costs$Chlorpromazine, 
            type = 'scatter', mode = 'markers', name = "Chlorpromazine",
            marker = list(symbol = "circle-open", color = 'black',size = 7)) %>%
  add_trace(x = ~psa_res$outcomes$Risperidone, y = ~psa_res$costs$Risperidone,
            type = 'scatter', mode = 'markers', name = "Risperidone",
            marker = list(symbol = "diamond-open", color = 'blue',size = 8)) %>%
  add_trace(x = ~psa_res$outcomes$Haloperidol, y = ~psa_res$costs$Haloperidol,
            type = 'scatter', mode = 'markers', name = "Haloperidol",
            marker = list(symbol = "x-thin-open", color = 'red',size = 8)) %>%
  add_trace(x = ~psa_res$outcomes$Quetiapine, y = ~psa_res$costs$Quetiapine,
            type = 'scatter', mode = 'markers', name = "Quetiapine",
            marker = list(symbol = "asterisk-open", color = 'orange',size = 8)) %>%
  add_trace(x = ~psa_res$outcomes$Olanzapine, y = ~psa_res$costs$Olanzapine,
            type = 'scatter', mode = 'markers', name = "Olanzapine",
            marker = list(symbol = "triangle-up-open", color = 'maroon',size = 8)) %>%
  layout(xaxis = list(title = 'DALYs averted'), yaxis = list(title = 'Costs'))



## COST EFFECTIVENESS ACCEPTBILITY CURVE
ceac_res <- RunCEAC(psa_output = psa_res, wtp = 4000, by = 1000)

plot_ly(ceac_res, x = ~WTP) %>%
  add_trace(y = ~Chlorpromazine, type = 'scatter', mode = 'lines+markers', name = "Chlorpromazine",
            line = list(shape = 'spline', color = col[1], width= 1.5),
            marker = list(symbol = "diamond", color = col[1], size = 9)) %>%
  
  add_trace(y = ~Haloperidol, type = 'scatter', mode = 'lines+markers', name = "Haloperidol",
            line = list(shape = 'spline', color = col[2], width= 1.5),
            marker = list(symbol = "square", color = col[2], size = 8)) %>%
  
  add_trace(y = ~Quetiapine, type = 'scatter', mode = 'lines+markers', name = "Quetiapine",
            line = list(shape = 'spline', color = col[3], width= 1.5),
            marker = list(symbol = "circle", color = col[3], size = 8)) %>%
  
  add_trace(y = ~Olanzapine, type = 'scatter', mode = 'lines+markers', name = "Olanzapine",
            line = list(shape = 'spline', color = col[4], width= 1.5),
            marker = list(symbol = "x-thin-open", color = col[4], size = 8)) %>%

  add_trace(y = ~Risperidone, type = 'scatter', mode = 'lines+markers', name = "Risperidone",
            line = list(shape = 'spline', color = col[5], width= 1.5),
            marker = list(symbol = "asterisk-open", color = col[5], size = 8)) %>%
  
  layout(xaxis = list(title = 'Willingness-to-pay ($/DALY averted)'), 
         yaxis = list(title = 'Probability Cost-Effective'),
         legend = list(orientation = 'v'))



# 
# psa_res <- lapply(split(psa_input, 1:nrow(psa_input)), FUN = RunPSA, base_case = model) %>%
#   do.call("rbind", .)
# 


# psa_res2 <- psa_input %>% 
#   purrr::transpose() %>% ## convert each row into a dataframe and combine dataframes into a list
#   map(.f = RunPSA, base_case = model) %>% # apply function to each dataframe in the list
#   bind_rows()



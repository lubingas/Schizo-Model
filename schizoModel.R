rm(list=ls())

# place all files in the same working directory and change to this directory 
setwd("../Schizo-Model/")

# load functions we need
source("sub_routines.R")

# some preliminary stuff
# these are not varied in sensitivity analyses
nStates <- 10
nCycles <- 70
interventions <- c("Chlorpromazine", "Haloperidol", "Quetiapine", "Risperidone", "Olanzapine")

source("data_input.R")

# pull out inputs for the main model but separate out the transition probabilities
baseinputs <- data %>%
  filter(!is.na(input) & dist != "dirichlet") %>%   # remove NA from data
  pull(input, variable)

## apply inflation adjustment where necessary

# process dirichlet inputs so I can input set of high and low values at the same time (also add label, so it is easy to reference in functions)
basetransitions <- data %>%
  filter(!is.na(input) & dist == "dirichlet") %>%   # remove NA from data
  pull(input, variable)

## run model
base_case <- RunModel(intervention = interventions, 
                      baseparms = baseinputs, 
                      basetransitions = basetransitions,
                      DALYs = TRUE)

# perform incremental analysis (establish efficiency frontier)
CEAnalysis <- IncrementalAnalysis(model = base_case)

# plot incremental analysis
ScatterPlot(CEAnalysis)

## ONE WAY SENSITIVITY ANALYSES
# parameters
usa.low <- data %>%
  filter(!is.na(low) & dist != "dirichlet") %>%   # remove NA from data
  pull(low, variable)

usa.high <- data %>%
  filter(!is.na(low) & dist != "dirichlet") %>%   # remove NA from data
  pull(high, variable)

# transition probabilities
dir_lowinputs <- data %>%
  filter(dist == "dirichlet") %>%   # remove NA from data
  pull(low, variable) %>%
  split(.,cut(seq_along(.), (length(.)/length(interventions)), labels = interventions))

dir_highinputs <- data %>%
  filter(dist == "dirichlet") %>%   # remove NA from data
  pull(high, variable) %>%
  split(.,cut(seq_along(.), (length(.)/length(interventions)), labels = interventions))

# iterate over variables and perform OWSA
dsa_res <- owsa(model = base_case, 
                low_base = usa.low, low_transitions = dir_lowinputs,
                high_base = usa.high, high_transitions = dir_highinputs,
                outcome = "cost", comparator = "Haloperidol")

# plot OWSA
TornadoPlot(owsa_list = dsa_res, max_vars = 10)


## PSA
# how many simulations?
nsims <- 10000


# first identify fixed parameters
fixedparms <- data %>%
  filter(dist == "fixed") %>%
  pull(variable)

betaparms <- data %>%
  filter(dist == "beta") %>%
  pull(variable)

normparms <- data %>%
  filter(dist == "normal") %>%
  pull(variable)

lognormparms <- data %>%
  filter(dist == "lognormal") %>%
  pull(variable)

alphas <- data %>%
  pull(alpha, variable)

betas <- data %>%
  pull(beta, variable)

sems <- data %>%
  pull(sem, variable)

# draw values
#betas
randbeta <- mapply(rbeta, MoreArgs=list(n=nsims), alphas[betaparms], betas[betaparms])

# normals
randnormal <- mapply(rnorm, MoreArgs=list(n=nsims), baseinputs[normparms], sems[normparms])

# log-normals
randlognormal <- mapply(rlnorm, MoreArgs=list(n=nsims), 
                           log(baseinputs[lognormparms]), 
                           (log(usa.high[lognormparms]) - log(usa.low[lognormparms]))/3.92)

# random dirichlet transitions
randdirchinputs <- data %>%
  filter(dist == "dirichlet") %>%   # remove NA from data
  pull(alpha, variable) %>%
  split(.,cut(seq_along(.), (length(.)/length(interventions)), labels = interventions))

# generate variates and put in one database
randdirichlets <- lapply(randdirchinputs, vdirichlet, n=nsims) %>%
  do.call("cbind", .)

fixedinputs <- matrix(rep(baseinputs[fixedparms],each=nsims),nrow=nsims) %>%
  magrittr::set_colnames(., fixedparms)

# pool together into one dataframe
psa_input <- data.frame(cbind(randbeta, randnormal, randlognormal, randdirichlets, fixedinputs))


psa_res <- RunPSA(base_case = base_case, psalist = psa_input)

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
ceac_res <- RunCEAC(costs = psa_res$costs, outcomes = psa_res$outcomes, wtp = 2000, by = 500)

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



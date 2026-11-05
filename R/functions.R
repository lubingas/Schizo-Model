# R/functions.R
# ---------------------------
# Core functions for the Schizophrenia cost-effectiveness model
# ---------------------------

## functions for calculating various stuff
# 1) extrapyramidal side effects
eps.risk <- function(base.p, drug.OR){
  return(drug.OR * base.p/(1 - base.p + (drug.OR * base.p)))
}

# 2) Diabetes risk
diab.risk <- function(base.p, drug.RR) {
  return(base.p * drug.RR)
}

# 3) Price inflation: ratio of the CPI in year `to` to the CPI in year `from`,
#    from the CPI sheet of data_input.xlsx (inputs$cpi). Costs are reported in 2012 US$.
cpi_ratio <- function(cpi, from, to = 2012) {
  cpi$Value[cpi$Year == to] / cpi$Value[cpi$Year == from]
}

# 4) costs models
# a) drug costs (all drug prices in data_input.xlsx are 2012 US$)
# inpatient: daily dose for the length of stay
ip.drug.costs <- function(cost, qty, days) {
  return(cost * qty * days)
}

op.drug.costs <- function(cost, qty) {
  return(cost * qty * 365)
}

## extract parameters from data
getparms <- function(data, value) {
  data %>% 
    filter(!is.na({{value}})) %>%   # remove NA from data
    select(variable, {{value}}) %>%
    deframe() # return named vector
}

# Option A (cleanest): infer drug from variable suffix
getTPset <- function(data, value) {
  data %>%
    filter(!is.na({{value}})) %>%
    select(variable, {{value}}) %>%
    deframe() %>%
    split(stringr::str_extract(names(.), "(CPZ|HAL|OLA|RIS|QUE)$"))
}

## function to get the low or high value of each transition probability for the one-way analysis
## INSTRUCTIONS
# 1) data: the transition probability table (inputs$transitions)
# 2) q: 0.025 for the low values, 0.975 for the high values
# each probability's marginal distribution is Beta(a_i, sum(a) - a_i), where a = the counts (alpha)
# of that drug's set, as in the Dirichlet draws of the PSA; returns a list by drug, like getTPset()
getTPbounds <- function(data, q) {
  data %>%
    filter(dist == "dirichlet") %>%
    getTPset(value = alpha) %>%
    lapply(function(a) qbeta(q, a, sum(a) - a))
}

## ----------------------------
## Pieces of the model that don't depend on the drug sequence
## Each takes the parameters it needs as arguments: params is the named vector of model
## parameters (base case, one-way or PSA draw), read by name with [[ ]]
## ----------------------------

## function to look up a table's value for the age of the cohort in each cycle
## INSTRUCTIONS
# 1) tab: a table by age band, with columns ages (first age of the band) and index (band number)
# 2) values: the value of each band to return, e.g. tab$p.death
# 3) ages: the cohort's age in each cycle
by_age <- function(tab, values, ages) {
  values[match(findInterval(x = ages, vec = tab$ages), tab$index)]
}

## function to collect one drug's parameters in the order mtrace() expects:
## c(disAE, disEff, disOther, pHosp, pWt, eps, diab, cost, qty)
## INSTRUCTIONS
# 1) params: named vector of model parameters
# 2) code: the drug's code at the end of its parameter names, e.g. "CPZ"
# 3) name: the drug's name at the start of its price and dose names, e.g. "chlorpromazine"
drug_parms <- function(params, code, name) {
  unname(params[c(paste0(c("disAE", "disEff", "disOther", "pHosp", "pWt", "eps", "diab"), code),
                  paste0(name, c("_cost", "_qty")))])
}

## function to work out the indirect cost of a carer's lost work time, per cycle (in hospital and as an outpatient)
## INSTRUCTIONS
# 1) params: named vector of model parameters
carer_indirect <- function(params) {
  params[["laborforcerate"]]*(1-params[["unemploymentrate"]])*params[["GDP_percapita"]]*params[["time"]]*params[["DALY_caretaker"]]*(264/365)
}

## function to work out the costs of a cycle in hospital that don't depend on the drug
## INSTRUCTIONS
# 1) params: named vector of model parameters
# 2) cpi: the CPI table (inputs$cpi), to bring 2001 lab prices to 2012
# returns a named list: personnel, lab, capital, carer.travel, patient.travel, indirect, carer.indirect
ip_costs <- function(params, cpi) {
  
  # b) personnel: staff time on ward rounds (per week, so divided by 7 for a day) and giving drugs (per day),
  #    over the length of stay. Each term is rounds x staff x hours x hourly pay (salary / work hours in a year).
  #    (The Excel model divided only the psychiatrist's major round by 7, and costed the drug round in hours, not US$)
  psy.major  <- params[["major_week"]]*params[["PSY_major"]]*(params[["duration_major"]]/60)*(params[["rate_psy"]]/params[["workhours_year"]])
  pco.major  <- params[["major_week"]]*params[["PCO_major"]]*(params[["duration_major"]]/60)*(params[["rate_pco"]]/params[["workhours_year"]])
  pno.major  <- params[["major_week"]]*params[["PNO_major"]]*(params[["duration_major"]]/60)*(params[["rate_pno"]]/params[["workhours_year"]])
  pco.minor  <- params[["minor_week"]]*params[["PCO_minor"]]*(params[["duration_minor"]]/60)*(params[["rate_pco"]]/params[["workhours_year"]])
  pno.minor  <- params[["minor_week"]]*params[["PNO_minor"]]*(params[["duration_minor"]]/60)*(params[["rate_pno"]]/params[["workhours_year"]])
  drug.admin <- params[["PNO_drugadmin"]]*(params[["duration_drugadmin"]]/60)*(params[["rate_pno"]]/params[["workhours_year"]])
  personnel  <- ((psy.major + pco.major + pno.major + pco.minor + pno.minor)/7 + drug.admin)*params[["time"]]*params[["LOS"]]
  
  # c) laboratory tests (2001 prices), done once per hospital stay: blood counts, syphilis (TPHA or VDRL) and HIV
  #    (the Excel model costed them 12 times a cycle, as if monthly for a year)
  tests <- params[["FBC_cost"]]*params[["FBC_freq"]] + params[["ESR_cost"]]*params[["ESR_freq"]] +
    params[["TPHA_cost"]]*params[["TPHA_freq"]]*params[["TPHA_proportion"]] +
    params[["VDRL_cost"]]*params[["VDRL_freq"]]*(1-params[["TPHA_proportion"]]) + params[["HIV_cost"]]*params[["HIV_freq"]]
  lab <- cpi_ratio(cpi, 2001)*tests
  
  # d) capital and overheads, per day in hospital
  capital <- params[["capital_costs_per_inpatient_per_day"]]*params[["time"]]*params[["LOS"]]
  
  # e) travel: carers' visits during the stay, and the patient's trip there and back
  carer.travel   <- params[["averagetravel_cost"]]*2*params[["caregiver_prop"]]*(params[["time"]]*params[["LOS"]]*params[["caregiver_freq"]])/params[["every"]]
  patient.travel <- params[["averagetravel_cost"]]*2
  
  # f) indirect: the patient's lost work time
  indirect <- params[["laborforcerate"]]*(1-params[["unemploymentrate"]])*params[["GDP_percapita"]]*params[["time"]]
  
  list(personnel = personnel, lab = lab, capital = capital, carer.travel = carer.travel,
       patient.travel = patient.travel, indirect = indirect, carer.indirect = carer_indirect(params))
}

## function to work out the costs of a cycle as an outpatient that don't depend on the drug
## INSTRUCTIONS
# 1) params: named vector of model parameters
# 2) cpi: the CPI table (inputs$cpi), to bring 2001 lab prices to 2012
# returns a named list: personnel, lab, capital, carer.travel, patient.travel, indirect, carer.indirect
op_costs <- function(params, cpi) {
  
  # b) personnel: the consultation (psychiatrist or clinical officer, and a nurse) and dispensing, per visit.
  #    Each term is staff x share of patients they see x hourly pay; the consultation length applies to all three
  #    (the Excel model applied it to the psychiatrist only, costing the others for a full hour)
  psy        <- params[["PSY"]]*params[["PSY_prop"]]*(params[["rate_psy"]]/params[["workhours_year"]])
  pco        <- params[["PCO"]]*(1-params[["PSY_prop"]])*(params[["rate_pco"]]/params[["workhours_year"]])
  nurse      <- params[["nurse"]]*params[["nurse_prop"]]*(params[["rate_pno"]]/params[["workhours_year"]])
  pharmacist <- params[["pharm"]]*(params[["duration_drugdisp"]]/60)*(params[["rate_pharm"]]/params[["workhours_year"]])
  personnel  <- ((params[["duration_consult"]]/60)*(psy + pco + nurse) + pharmacist)*params[["visitation_frequency"]]*params[["time"]]*12
  
  # c) laboratory tests (2001 prices): blood counts once a cycle
  #    (the Excel model costed them 12 times a cycle, as if monthly)
  lab <- cpi_ratio(cpi, 2001)*(params[["FBC_cost"]]*params[["FBC_freq"]] + params[["ESR_cost"]]*params[["ESR_freq"]])
  
  # d) capital and overheads, per visit
  capital <- params[["capital_costs_per_outpatient_visit"]]*params[["visitation_frequency"]]*params[["time"]]*12
  
  # e) travel to each visit: the carer (for the share of outpatients who come with one) and the patient
  #    (the Excel model used the share of hospital patients visited by a carer, caregiver_prop)
  carer.travel   <- params[["caregiver_visit"]]*params[["averagetravel_cost"]]*2*params[["visitation_frequency"]]*params[["time"]]*12
  patient.travel <- params[["averagetravel_cost"]]*2*params[["visitation_frequency"]]*params[["time"]]*12
  
  # f) indirect: the patient's lost work days for visits
  indirect <- params[["laborforcerate"]]*(1-params[["unemploymentrate"]])*(params[["GDP_percapita"]]/264)*(params[["visitation_frequency"]]*params[["time"]]*12)
  
  list(personnel = personnel, lab = lab, capital = capital, carer.travel = carer.travel,
       patient.travel = patient.travel, indirect = indirect, carer.indirect = carer_indirect(params))
}

## function to work out the disability weight of a cycle on a drug: the phase, plus weight gain, EPS and diabetes
## INSTRUCTIONS
# 1) params: named vector of model parameters
# 2) drug: the drug's parameters from drug_parms() (uses pWt, eps and diab)
# 3) disutilDIAB: disability weight of diabetes (worked out in RunModel)
# 4) phase: disability weight of the phase -- residualDALY, or acuteDALY in hospital
state_weight <- function(params, drug, disutilDIAB, phase = params[["residualDALY"]]) {
  phase + params[["disutilWeight"]]*drug[5] +
    eps.risk(params[["pEPSPLA"]], drug[6])*params[["durEPS"]]*params[["disutilEPS"]] +
    diab.risk(params[["pDiabetesHAL"]], drug[7])*disutilDIAB
}

## function to work out the disability weight of a cycle on fluphenazine decanoate (the 3rd-line drug)
## INSTRUCTIONS
# 1) params: named vector of model parameters
# 2) disutilDIAB: disability weight of diabetes (worked out in RunModel)
# fluphenazine has its own diabetes probability rather than a relative risk, so it doesn't use state_weight()
third_line_weight <- function(params, disutilDIAB) {
  params[["residualDALY"]] + params[["disutilWeight"]]*params[["pWtFLU"]] +
    eps.risk(params[["pEPSPLA"]], params[["epsFLU"]])*params[["durEPS"]]*params[["disutilEPS"]] +
    params[["pDiabetesFLU"]]*disutilDIAB
}

## function to work out the cost of a cycle as an outpatient on a drug: the drug, benzhexol for EPS,
## diabetes care, and the outpatient costs that don't depend on the drug
## INSTRUCTIONS
# 1) params: named vector of model parameters
# 2) drug: the drug's parameters from drug_parms() (uses eps, diab, cost and qty)
# 3) op: outpatient costs from op_costs()
residual_cost <- function(params, drug, op) {
  op.drug.costs(drug[8], drug[9]) +
    eps.risk(params[["pEPSPLA"]], drug[6])*op.drug.costs(params[["artane_cost"]], params[["artane_qty"]])*params[["durEPS"]] +
    diab.risk(params[["pDiabetesHAL"]], drug[7])*params[["cDiabetes"]] +
    op$personnel + op$lab + op$capital + op$carer.travel + op$patient.travel + op$indirect + op$carer.indirect
}

## function to work out the cost of a cycle as an outpatient on fluphenazine decanoate (the 3rd-line drug)
## INSTRUCTIONS
# 1) params: named vector of model parameters
# 2) op: outpatient costs from op_costs()
# a monthly injection, so priced per month; its own diabetes probability rather than a relative risk
third_line_cost <- function(params, op) {
  params[["fluphenazine_cost"]]*params[["fluphenazine_qty"]]*12 +
    eps.risk(params[["pEPSPLA"]], params[["epsFLU"]])*op.drug.costs(params[["artane_cost"]], params[["artane_qty"]])*params[["durEPS"]] +
    params[["pDiabetesFLU"]]*params[["cDiabetes"]] +
    op$personnel + op$lab + op$capital + op$carer.travel + op$patient.travel + op$indirect + op$carer.indirect
}

## function to work out the cost of a cycle with an acute episode: a hospital stay on a drug
## (chlorpromazine in the model), benzhexol for EPS, diabetes care, and the hospital costs that don't depend on the drug
## INSTRUCTIONS
# 1) params: named vector of model parameters
# 2) drug: the drug's parameters from drug_parms() (uses eps, diab, cost and qty)
# 3) ip: hospital costs from ip_costs()
acute_cost <- function(params, drug, ip) {
  ip.drug.costs(drug[8], drug[9], params[["LOS"]]) +
    eps.risk(params[["pEPSPLA"]], drug[6])*ip.drug.costs(params[["artane_cost"]], params[["artane_qty"]], params[["LOS"]]) +
    diab.risk(params[["pDiabetesHAL"]], drug[7])*params[["cDiabetes"]] +
    ip$personnel + ip$lab + ip$capital + ip$carer.travel + ip$patient.travel + ip$indirect + ip$carer.indirect
}


## function to run the Markov trace, DALYs and costs for one treatment sequence
## INSTRUCTIONS
# 1) primary, secondary: the first- and second-line drug's parameters, from drug_parms()
# 2) pRemit: probability of remission when treated in hospital
# 3) params: named vector of model parameters
# 4) shared: what RunModel() works out once for all five sequences -- the death probabilities, life
#    expectancy and discount factors for each cycle, and the disability weights and costs of the states
#    that don't depend on the drugs
# states: 1 residual on 1st line, 2 residual on 2nd line, 3 residual on 3rd line, 4 residual off 1st line,
#   5 residual off 2nd line, 6 acute off 1st line, 7 acute off 2nd line, 8 acute was on 1st line,
#   9 acute was on 2nd line, 10 dead
mtrace <- function(primary, secondary, pRemit, params, shared) {
  
  nStates <- shared$nStates
  nCycles <- shared$nCycles
  pDie <- shared$pDie
  
  ## transition probability array; the 3rd dimension holds the age-specific death probabilities
  trans <- array(0,c(nStates,nStates,nCycles))
  trans[1,2,] <- primary[1] + primary[2]
  trans[2,3,] <- secondary[1] + secondary[2]
  
  # residual on 1st line to residual off 1st line and residual on 2nd line to residual off 2nd line
  trans[1,4,] <- primary[3]
  trans[2,5,] <- secondary[3]
  
  # residual on 1st line to acute was on 1st line and residual on 2nd line to acute was on 2nd line
  trans[1,8,] <- primary[4]
  trans[2,9,] <- secondary[4]
  
  trans[4,6,] <- trans[5,7,] <- params[["pHospUntrt"]]
  
  trans[6,1,] <- trans[7,2,] <- trans[8,2,] <- trans[9,3,] <- pRemit
  
  ## get time varying death probabilities, for every cycle at once
  # insert death probabilities: column nStates of each cycle's matrix
  trans[,nStates,] <- rep(pDie, each = nStates)
  # row sums of every cycle's matrix, as a [state, cycle] matrix (the diagonal is still 0)
  rowsums <- colSums(aperm(trans, c(2, 1, 3)))
  # insert probability of remaining in the same state in each cycle (complement of rowsums)
  # idx picks the diagonal of every cycle's matrix: [1,1,1], [2,2,1], ..., [nStates,nStates,nCycles]
  idx <- cbind(rep(1:nStates, nCycles), rep(1:nStates, nCycles), rep(1:nCycles, each = nStates))
  trans[idx] <- 1 - rowsums
  # set absorbing state (death) probability to 1
  trans[nStates,nStates,] <- 1
  
  # cohort trace [cycle, state]: everyone starts in state 1 (residual, on the first-line drug)
  # and moves with the transition matrix for their age in that cycle
  pop.trace <- matrix(0, nCycles, nStates)
  pop.trace[1, ] <- c(1,0,0,0,0,0,0,0,0,0)
  for (c in 2:nCycles) {
    pop.trace[c, ] <- pop.trace[c-1, ] %*% trans[,,c-1]
  }
  
  # calculate YLL: in cycle i, those alive in cycle i-1 times pDie[i-1] times expLife[i] (none in cycle 1)
  int.YLL <- rowSums(pop.trace[,1:9])
  YLL <- c(0, int.YLL[-nCycles] * shared$expLife[-1] * pDie[-nCycles])
  
  # disability weight of each state: the two drugs of this sequence, then the states every sequence shares
  utils.vec <- c(state_weight(params, primary, shared$disutilDIAB),
                 state_weight(params, secondary, shared$disutilDIAB),
                 shared$weights.other)
  YLD <- pop.trace %*% utils.vec
  DALYS <- shared$discFactOutcomes * (YLL + YLD)
  
  # cost of a cycle in each state, in the same order
  cost.vec <- c(residual_cost(params, primary, shared$op),
                residual_cost(params, secondary, shared$op),
                shared$costs.other)
  costs <- shared$discFactCosts * (pop.trace %*% cost.vec)
  
  # trace: proportion of the cohort in each state in each cycle
  list(DALYS = sum(DALYS), Costs = sum(costs), trace = pop.trace)
}


## function to run the model: costs and DALYs of the five treatment sequences
## INSTRUCTIONS
# 1) inputs: from read_inputs() -- tables (life tables, CPI), nStates, nCycles and interventions
# 2) baseparms, basetransitions: named vectors of parameters and transition probabilities
#    (base case, one-way or PSA values)
# 3) DALYs: TRUE when outcomes are DALYs (fewer is better)
# 4) summary_only: TRUE returns just the costs and outcomes by drug (all the one-way and PSA runs need)
RunModel <- function(inputs,
                     baseparms,
                     basetransitions,
                     DALYs = TRUE,
                     summary_only = FALSE) {
  
  # model parameters: one named vector, read by name with [[ ]]
  params <- c(baseparms, basetransitions)
  
  # Background mortality: the life table gives the probability of dying within each age band
  # (5-year bands in the WHO life table for Uganda), so convert it to a yearly probability
  lt <- inputs$lifeTab
  band <- c(diff(lt$ages), tail(diff(lt$ages), 1))
  q_year <- 1 - (1 - lt$p.death)^(1 / band)
  
  ## everything below depends only on the parameters, so it is worked out once and shared by the five sequences
  cycles <- seq(0, inputs$nCycles-1, by = params[["time"]])
  
  # discount factor
  discFactCosts <- 1 / (1 + params[["discCosts"]]) ^ cycles
  discFactOutcomes <- 1 / (1 + params[["discOutcomes"]]) ^ cycles
  
  ## age in each cycle, matched with the probability of death and the life expectancy for that age
  ages <- params[["age"]] + cycles
  p.death <- by_age(inputs$lifeTab, q_year, ages)
  p.Schizo.death <- by_age(inputs$schizoTab, inputs$schizoTab$p.death, ages)
  expLife <- by_age(inputs$lifeExp, inputs$lifeExp$yrs, ages)
  
  # calculate probability of death (from other causes, schizophrenia specific -- but not both)
  pDie <- p.death + p.Schizo.death - (p.death * p.Schizo.death)
  
  # diabetes mean disutility (in RunModel since its inputs vary in the sensitivity analyses)
  disutilDIAB <- params[["pNeuro"]]*params[["dNeuro"]] + params[["pFoot"]]*params[["dFoot"]] +
    params[["pImpotent"]]*params[["dImpotent"]] + params[["pNephro"]]*params[["dNephro"]] + params[["dDiabetes"]]
  
  # costs per cycle that don't depend on the drug: in hospital and as an outpatient
  ip <- ip_costs(params, inputs$cpi)
  op <- op_costs(params, inputs$cpi)
  
  # the states every sequence shares (3 to 10), in mtrace()'s order: residual on the 3rd-line drug,
  # residual off the 1st- and 2nd-line drug, the four acute states (in hospital, on chlorpromazine), dead
  cpz <- drug_parms(params, "CPZ", "chlorpromazine")
  acute.weight <- state_weight(params, cpz, disutilDIAB, phase = params[["acuteDALY"]])
  acute.cost <- acute_cost(params, cpz, ip)
  
  shared <- list(nStates = inputs$nStates, nCycles = inputs$nCycles, pDie = pDie, expLife = expLife,
                 discFactCosts = discFactCosts, discFactOutcomes = discFactOutcomes,
                 disutilDIAB = disutilDIAB, op = op,
                 weights.other = c(third_line_weight(params, disutilDIAB),
                                   params[["residualDALY"]], params[["residualDALY"]],
                                   acute.weight, acute.weight, acute.weight, acute.weight, 0),
                 costs.other = c(third_line_cost(params, op), 0, 0,
                                 acute.cost, acute.cost, acute.cost, acute.cost, 0))
  
  ## Markov traces for each treatment sequence (in the order of inputs$interventions):
  ## first-line drug, second-line drug (fluphenazine decanoate is third line for all)
  hal <- drug_parms(params, "HAL", "haloperidol")
  que <- drug_parms(params, "QUE", "quetiapine")
  ris <- drug_parms(params, "RIS", "risperidone")
  ola <- drug_parms(params, "OLA", "olanzapine")
  traces <- list(
    mtrace(primary = cpz, secondary = ris, pRemit = params[["pRemitCPZ"]], params = params, shared = shared),
    mtrace(primary = hal, secondary = ris, pRemit = params[["pRemitCPZ"]], params = params, shared = shared),
    mtrace(primary = que, secondary = hal, pRemit = params[["pRemitCPZ"]], params = params, shared = shared),
    mtrace(primary = ris, secondary = hal, pRemit = params[["pRemitCPZ"]], params = params, shared = shared),
    mtrace(primary = ola, secondary = hal, pRemit = params[["pRemitCPZ"]], params = params, shared = shared)
  )
  
  outcomes <- sapply(traces, `[[`, "DALYS")
  costs <- sapply(traces, `[[`, "Costs")

  # costs and outcomes only, named by drug (all the one-way and probabilistic analyses need)
  if(summary_only) return(list(costs = setNames(costs, inputs$interventions),
                               outcomes = setNames(outcomes, inputs$interventions)))

  output_tab <- tibble(intervention = inputs$interventions,
                       costs = costs,
                       outcomes = outcomes)

  # return things in output I need for other functions
  return(list(ce_result = output_tab,
              baseparms = baseparms, 
              basetransitions = basetransitions, 
              DALYs = DALYs,
              inputs = inputs,
              traces = setNames(lapply(traces, `[[`, "trace"), inputs$interventions)))
  
}

# function to identify interventions on the efficiency frontier
# INSTRUCTIONS
# input 3 vectors: 1) intervention, 2) costs, and 3) outcomes
# if we are working with DALYs, then reverse sort order
IncrementalAnalysis <- function(model) {

  # Cost-effectiveness frontier (efficiency frontier):
  #   1. drop options that are strongly dominated (another option costs no more
  #      and is at least as effective, and is strictly better on one of the two)
  #   2. order the rest by cost and drop options that are extendedly dominated
  #      (their ICER is higher than that of the next, more effective option)
  #   3. report ICERs between consecutive options on the frontier
  res_tab <- model$ce_result
  DALYs   <- model$DALYs
  # effectiveness: higher is better (DALYs incurred are negated)
  eff  <- if (DALYs) -res_tab$outcomes else res_tab$outcomes
  cost <- res_tab$costs
  n    <- nrow(res_tab)

  strongly <- vapply(seq_len(n), function(i) {
    any(cost[-i] <= cost[i] & eff[-i] >= eff[i] & (cost[-i] < cost[i] | eff[-i] > eff[i]))
  }, logical(1))

  cand <- which(!strongly)
  cand <- cand[order(cost[cand], -eff[cand])]
  extended <- integer(0)
  repeat {
    if (length(cand) < 3) break
    icer <- diff(cost[cand]) / diff(eff[cand])
    bad  <- which(icer[-length(icer)] > icer[-1])
    if (!length(bad)) break
    extended <- c(extended, cand[bad[1] + 1])
    cand <- cand[-(bad[1] + 1)]
  }
  frontier <- res_tab$intervention[cand]

  status <- ifelse(strongly, "Dominated",
            ifelse(seq_len(n) %in% extended, "Extendedly dominated", "On frontier"))

  # all options, relative to the least effective one (used for the scatter plot)
  ce_points <- res_tab %>%
    mutate(status = status) %>%
    arrange(if (DALYs) desc(outcomes) else outcomes) %>%
    mutate(delta_outcomes = if (DALYs) first(outcomes) - outcomes else outcomes - first(outcomes),
           delta_costs = costs - first(costs))

  ce_lines <- ce_points %>%
    filter(intervention %in% frontier) %>%
    arrange(delta_outcomes)

  # ICERs between consecutive frontier options, cheapest first
  incrementals <- ce_lines %>%
    select(intervention, costs, outcomes) %>%
    arrange(costs) %>%
    mutate(delta_costs    = costs - lag(costs),
           delta_outcomes = if (DALYs) lag(outcomes) - outcomes else outcomes - lag(outcomes),
           iCER           = delta_costs / delta_outcomes,
           explanation    = if_else(is.na(lag(intervention)),
                                    if (n() == 1) "Ref (dominates all other options)" else "Ref",
                                    paste0("iCER versus ", lag(intervention)))) %>%
    mutate(across(c(delta_costs, delta_outcomes, iCER), ~ replace(.x, is.na(.x), 0)))

  list(ce_points = ce_points,
       ce_lines = ce_lines,
       incrementals = incrementals,
       dominated = res_tab$intervention[strongly],
       extended_dominated = res_tab$intervention[extended],
       outcome_to_plot = if (DALYs) "DALYs averted" else "QALYs gained")
}


ScatterPlot <- function(FrontierAnalysis) {
  
  plot_ly() %>%
    add_trace(data = FrontierAnalysis$ce_points, x = ~delta_outcomes, y = ~delta_costs, 
              alpha = 0.7,
              color = ~factor(intervention),
              stroke = ~factor(intervention),
              #symbol = ~factor(intervention),
              size = 20,
              type="scatter", mode="markers") %>%
    add_trace(data = FrontierAnalysis$ce_lines, x = ~delta_outcomes, y = ~delta_costs, 
              type = "scatter", mode="lines", name = "CEEF", 
              line = list(color = 'rgb(22, 96, 167)', width = 2, dash = "dot")) %>%
    layout(yaxis = list(title = "Incremental Costs (US$)",
                        ticks = "outside"),
           xaxis = list(title = FrontierAnalysis$outcome_to_plot,
                        ticks = "outside"),
           font = list(size=12, 
                       family='Fira Sans'))
  
}


## function to set one probability of a drug's transition set and rescale the others proportionally,
## so the set still sums to 1
## INSTRUCTIONS
# 1) shares: named vector of the base-case probabilities in one set (sums to 1)
# 2) n: name of the probability to change
# 3) value: its new value
rebalance <- function(shares, n, value) {
  others <- names(shares) != n
  shares[others] <- shares[others] * (1 - value) / (1 - shares[n])
  shares[n] <- value
  return(shares)
}


## function to run the one-way sensitivity analysis
## INSTRUCTIONS
# 1) model: the base-case model object from RunModel()
# 2) low_base, high_base: low and high values of the parameters, from getparms()
# 3) low_transitions, high_transitions: low and high value of each transition probability, by drug,
#    from getTPbounds()
# each run starts from a fresh copy of the base case with one variable changed, and keeps the costs
# and outcomes of every drug; owsaPair() then picks the drugs and the result to show
owsa <- function(model, 
                 low_base, low_transitions, 
                 high_base, high_transitions) {
  
  # get the base-case inputs from the model object
  baseparms <- model$baseparms
  basetransitions <- model$basetransitions
  base <- c(baseparms, basetransitions)
  
  # 1) task list: for each variable, the values that change in its low run and in its high run
  # a parameter changes only itself
  baseTasks <- lapply(names(low_base), function(n) {
    list(varname = n, low = low_base[n], high = high_base[n])
  })
  
  # a transition probability changes its drug's whole set: itself, and the others rescaled
  transitionTasks <- lapply(names(low_transitions), function(s) {
    shares <- basetransitions[names(low_transitions[[s]])]   # base-case probabilities for this drug
    lapply(names(shares), function(n) {
      list(varname = n,
           low  = rebalance(shares, n, low_transitions[[s]][n]),
           high = rebalance(shares, n, high_transitions[[s]][n]))
    })
  })
  
  # one list of tasks: parameters, then every transition probability of every drug
  tasks <- c(baseTasks, unlist(transitionTasks, recursive = FALSE))
  
  # 2) costs and outcomes of every drug from one run
  runLive <- function(live) {
    RunModel(model$inputs,
             baseparms = live[names(baseparms)],
             basetransitions = live[names(basetransitions)],
             DALYs = model$DALYs, summary_only = TRUE)
  }
  
  # 3) run each task: a fresh copy of the base case for its low run and for its high run
  results <- lapply(seq_along(tasks), function(i) {
    task <- tasks[[i]]
    message("Running DSA on ", task$varname, ", variable ", i, " of ", length(tasks))
    # a misspelled name would otherwise be added to live silently
    stopifnot(all(c(names(task$low), names(task$high)) %in% names(base)))
    
    live_low <- base
    live_low[names(task$low)] <- task$low
    
    live_high <- base
    live_high[names(task$high)] <- task$high
    
    list(low = runLive(live_low), high = runLive(live_high))
  })
  
  # 4) one matrix per result: a row per variable, a column per drug
  varnames <- sapply(tasks, function(task) task$varname)
  collect <- function(run, what) {
    out <- do.call(rbind, lapply(results, function(r) r[[run]][[what]]))
    rownames(out) <- varnames
    return(out)
  }
  
  basecase <- runLive(base)
  
  return(list(variable = varnames,
              low_in = sapply(tasks, function(task) task$low[[task$varname]]),
              high_in = sapply(tasks, function(task) task$high[[task$varname]]),
              costs_low = collect("low", "costs"), costs_high = collect("high", "costs"),
              outcomes_low = collect("low", "outcomes"), outcomes_high = collect("high", "outcomes"),
              costs_base = basecase$costs, outcomes_base = basecase$outcomes,
              DALYs = model$DALYs))
  
}


## function to pick one result of the one-way analysis: intervention versus comparator
## INSTRUCTIONS
# 1) owsa_res: output of owsa()
# 2) outcome: "cost" (incremental cost), "outcome" (DALYs averted / QALYs gained), "ICER",
#    or "NMB" (incremental net monetary benefit: outcome valued at wtp, minus incremental cost)
# 3) intervention, comparator: two of the drugs
# 4) wtp: willingness to pay per DALY averted / QALY gained; needed for "NMB"
# returns the table TornadoPlot() uses; no model runs, so it is instant
owsaPair <- function(owsa_res, outcome, intervention, comparator, wtp = NULL) {
  
  stopifnot(outcome %in% c("cost", "outcome", "ICER", "NMB"), outcome != "NMB" || is.numeric(wtp))
  
  # Result of interest: `intervention` versus `comparator`;
  # costs and outcomes are matrices with a row per run and a column per drug
  pairwise <- function(costs, outcomes) {
    dc <- costs[, intervention] - costs[, comparator]
    do <- if (owsa_res$DALYs) outcomes[, comparator] - outcomes[, intervention]
          else outcomes[, intervention] - outcomes[, comparator]
    value <- unname(switch(outcome, ICER = dc / do, cost = dc, outcome = do, NMB = wtp * do - dc))
    # an ICER is a cost per unit of health only when the intervention costs more and gains more;
    # in any other quadrant the ratio can't be ranked, so say which quadrant instead
    status <- rep(NA_character_, length(value))
    if (outcome == "ICER") {
      status <- unname(ifelse(dc > 0 & do > 0, NA,
                       ifelse(dc <= 0 & do >= 0, "dominant",
                       ifelse(dc >= 0 & do <= 0, "dominated", "cheaper, less effective"))))
    }
    list(value = value, status = status)
  }
  
  base.case <- pairwise(rbind(owsa_res$costs_base), rbind(owsa_res$outcomes_base))
  low.run <- pairwise(owsa_res$costs_low, owsa_res$outcomes_low)
  high.run <- pairwise(owsa_res$costs_high, owsa_res$outcomes_high)
  
  owsavars <- paste0(owsa_res$variable, " (", sprintf(owsa_res$low_in, fmt = "%#.2f"), ", ",
                     sprintf(owsa_res$high_in, fmt = "%#.2f"), ")")
  
  # post processing of data to return full results; for the ICER, runs that leave the
  # cost-more-gain-more quadrant have no ratio (NA) and sort above the rest, since they change the conclusion
  owsa_tab <- tibble(variable = owsavars,
                     low = ifelse(is.na(low.run$status), low.run$value, NA),
                     high = ifelse(is.na(high.run$status), high.run$value, NA),
                     low_status = low.run$status, high_status = high.run$status) %>%
    mutate(high_val = high - base.case$value, low_val = low - base.case$value,
           range = coalesce(abs(high - low), abs(low_val), abs(high_val), 0),
           no_ratio = (!is.na(low_status)) + (!is.na(high_status))) %>%
    arrange(no_ratio, range)
  
  return(list(owsa_data = owsa_tab, owsa_outcome = outcome, owsa_basecase = base.case$value,
              owsa_basestatus = base.case$status, wtp = wtp, DALYs = owsa_res$DALYs))
  
}

## function to say, for each row of a tornado, which runs ended up in the given quadrants (ICER only)
## e.g. "low: dominant", "high: dominated", "low and high: dominant"; NA when neither did
tornadoLabels <- function(owsa_tab, statuses) {
  lo <- ifelse(owsa_tab$low_status %in% statuses, owsa_tab$low_status, NA)
  hi <- ifelse(owsa_tab$high_status %in% statuses, owsa_tab$high_status, NA)
  case_when(!is.na(lo) & !is.na(hi) & lo == hi ~ paste("low and high:", lo),
            !is.na(lo) & !is.na(hi) ~ paste0("low: ", lo, "; high: ", hi),
            !is.na(lo) ~ paste("low:", lo),
            !is.na(hi) ~ paste("high:", hi),
            TRUE ~ NA_character_)
}

# function takes output from DSA as a list, and plots the results (user inputs max vars)
TornadoPlot <- function(owsa_list, max_vars) {
  
  # extract the things we need
  owsa_tab <- owsa_list$owsa_data
  outcome <- owsa_list$owsa_outcome
  basecase <- owsa_list$owsa_basecase
  DALYs <- owsa_list$DALYs
  on_outcome <- ifelse(DALYs == TRUE, "DALYs averted", "QALYs gained")
  if (outcome == "ICER" && !is.na(owsa_list$owsa_basestatus)) {
    stop("At the base case the intervention is ", owsa_list$owsa_basestatus,
         ", so there is no ICER to vary; use outcome = \"cost\", \"outcome\" or \"NMB\"")
  }
  
  # set color
  col_low <- "rgb(166,206,227)"
  col_high <- "rgb(31,120,180)"
  
  # owsa_tab comes sorted by effect, largest last; draw each row's longer bar first,
  # so the shorter one stays visible when both fall on the same side of the base case
  owsa_torn <- owsa_tab %>%
    slice_tail(n = max_vars) %>%
    mutate(low_behind = is.na(high_val) | (!is.na(low_val) & abs(low_val) >= abs(high_val)),
           back = as.numeric(ifelse(low_behind, low_val, high_val)), front = as.numeric(ifelse(low_behind, high_val, low_val)),
           back_col = ifelse(low_behind, col_low, col_high), front_col = ifelse(low_behind, col_high, col_low))
  
  axislabel <- switch(outcome,
                      ICER=paste0("Cost per ", on_outcome),
                      cost="Incremental Cost ($)",
                      outcome=on_outcome,
                      NMB=paste0("Net monetary benefit at $", format(owsa_list$wtp, big.mark = ","),
                                 " per ", sub("s .*", "", on_outcome), " ($)"))
  
  # ICER only: runs with no ratio are named at the edges, dominant on the left, the rest on the right
  xrange <- NULL; notes <- list()
  if (outcome == "ICER") {
    xr <- range(c(basecase, owsa_torn$low, owsa_torn$high), na.rm = TRUE)
    pad <- 0.45 * max(diff(xr), abs(basecase))
    xrange <- c(xr[1] - pad, xr[2] + pad)
    left <- tornadoLabels(owsa_torn, "dominant")
    right <- tornadoLabels(owsa_torn, c("dominated", "cheaper, less effective"))
    for (i in seq_len(nrow(owsa_torn))) {
      if (!is.na(left[i])) notes <- c(notes, list(list(x = xrange[1], y = owsa_torn$variable[i], text = left[i],
                                                        xanchor = "left", showarrow = FALSE)))
      if (!is.na(right[i])) notes <- c(notes, list(list(x = xrange[2], y = owsa_torn$variable[i], text = right[i],
                                                         xanchor = "right", showarrow = FALSE)))
    }
  }
  
  plot_ly(data = owsa_torn, y = ~variable) %>%
    add_trace(x = ~back, base = basecase, type = 'bar', orientation = 'h', showlegend = FALSE,
              marker = list(color = ~back_col)) %>%
    add_trace(x = ~front, base = basecase, type = 'bar', orientation = 'h', showlegend = FALSE,
              marker = list(color = ~front_col)) %>%
    # legend entries only (zero-length bars)
    add_trace(x = 0, y = owsa_torn$variable[1], base = basecase, type = 'bar', orientation = 'h',
              name = "lower input", marker = list(color = col_low), hoverinfo = "skip") %>%
    add_trace(x = 0, y = owsa_torn$variable[1], base = basecase, type = 'bar', orientation = 'h',
              name = "higher input", marker = list(color = col_high), hoverinfo = "skip") %>%
    layout(barmode='overlay') %>% 
    layout(yaxis = list(categoryorder = "array",
                        categoryarray = ~variable,
                        title = ""),
           xaxis = list(zeroline = FALSE,
                        showline = FALSE,
                        showgrid = FALSE,
                        range = xrange,
                        title = axislabel),
           annotations = notes,
           legend = list(itemclick = FALSE, itemdoubleclick = FALSE),
           font = list(size=12, 
                       family='Fira Sans'))
  
}


## ----------------------------
## PSA utilities
## ----------------------------

## function to draw every PSA input at once -- called once per analysis run
## INSTRUCTIONS
# 1) nsims: number of draws
# 2) parameters, transitions: the input tables (dist, alpha, beta, sem and the base case, low and high columns)
# 3) interventions: the drugs (kept for the call signature; the transition sets are found from the variable names)
# returns a numeric matrix: a row per draw, a named column per parameter
draw_psa_df <- function(nsims, parameters, transitions, interventions,
                        basecase = "basecase", low = "low", high = "high") {
  
  # --- ONE-TIME parameter classification i.e., parameters by distribution type ----
  betaparms    <- parameters %>% filter(dist == "beta")      %>% pull(variable)
  normparms    <- parameters %>% filter(dist == "normal")    %>% pull(variable)
  lognormparms <- parameters %>% filter(dist == "lognormal") %>% pull(variable)
  fixedparms   <- parameters %>% filter(dist == "fixed")     %>% pull(variable)
  
  alphas <- parameters %>% pull(alpha, name = variable)
  betas  <- parameters %>% pull(beta,  name = variable)
  sems   <- parameters %>% pull(sem,   name = variable)
  
  base_vals <- getparms(parameters, value = basecase)
  low_vals  <- getparms(parameters, value = low)
  high_vals <- getparms(parameters, value = high)
  
  # one column of nsims draws per parameter, in the order given
  draw <- function(parms, f) {
    matrix(vapply(parms, f, numeric(nsims)), nrow = nsims, dimnames = list(NULL, parms))
  }
  
  # --- Beta, Normal, Lognormal and fixed draws
  beta_m   <- draw(betaparms, function(v) rbeta(nsims, alphas[[v]], betas[[v]]))
  normal_m <- draw(normparms, function(v) rnorm(nsims, mean = base_vals[[v]], sd = sems[[v]]))
  logn_m   <- draw(lognormparms, function(v) rlnorm(nsims, meanlog = log(base_vals[[v]]),
                                                    sdlog = (log(high_vals[[v]]) - log(low_vals[[v]])) / 3.92))
  fixed_m  <- draw(fixedparms, function(v) rep(base_vals[[v]], nsims))
  
  # --- Dirichlet draws: one set per drug (drug code at the end of the name), in table order
  dir_alpha <- transitions %>% filter(dist == "dirichlet") %>% pull(alpha, variable)
  drug <- stringr::str_extract(names(dir_alpha), "(CPZ|HAL|OLA|RIS|QUE)$")
  dir_m <- do.call(cbind, lapply(split(dir_alpha, factor(drug, levels = unique(drug))), function(alpha) {
    k <- length(alpha)
    G <- matrix(rgamma(nsims * k, shape = rep(alpha, each = nsims)), nrow = nsims)
    D <- G / rowSums(G)
    colnames(D) <- names(alpha)
    D
  }))
  
  # --- Combine all blocks
  return(cbind(beta_m, normal_m, logn_m, dir_m, fixed_m))
  
}

## function to run the model on one PSA draw
## INSTRUCTIONS
# 1) model: the base-case model object from RunModel()
# 2) psalist: one draw, a named numeric vector of every parameter
genPSA <- function(model, psalist) {
  
  # get things we need
  basenames <- names(model$baseparms)
  transitnames <- names(model$basetransitions)
  stopifnot(all(c(basenames, transitnames) %in% names(psalist)))   # every parameter has a draw
  
  # costs and outcomes of every drug, named by drug
  res <- RunModel(model$inputs,
                  baseparms = psalist[basenames],
                  basetransitions = psalist[transitnames],
                  DALYs = model$DALYs, summary_only = TRUE)

  return(res)
  
}

## function to run the probabilistic sensitivity analysis
## INSTRUCTIONS
# 1) model: the base-case model object from RunModel()
# 2) nsims: number of draws
# 3) parameters, transitions: the tables to draw from (by default, the model's own inputs)
# progress is reported with progressr: wrap the call in with_progress() (console) or withProgressShiny() (app)
RunPSA <- function(model, nsims,
                   parameters = model$inputs$parameters,
                   transitions = model$inputs$transitions) {
  
  # draw every input once, then split the draws into rows: each row is a named numeric vector of every parameter
  psalist <- draw_psa_df(nsims = nsims, parameters = parameters, transitions = transitions,
                         interventions = model$inputs$interventions)
  draws <- asplit(psalist, MARGIN = 1)
  
  # update the progress bar 100 times, not once per draw
  every <- ceiling(length(draws) / 100)
  p <- progressr::progressor(steps = ceiling(length(draws) / every))
  
  results <- lapply(seq_along(draws), function(i) {
    if (i %% every == 0 || i == length(draws)) {
      p(message = paste0("PSA iteration ", i, " of ", length(draws)))
    }
    genPSA(model = model, psalist = draws[[i]])
  })
  
  # one table per result: a row per draw, a column per drug
  costs    <- as_tibble(do.call(rbind, unname(lapply(results, `[[`, "costs"))))
  outcomes <- as_tibble(do.call(rbind, unname(lapply(results, `[[`, "outcomes"))))
  
  return(list(costs = costs, outcomes = outcomes))
}

# function takes the maxwtp, a vector or dataframe of costs and outcomes for all interventions
# and returns probability that a given intervention is cost effective
# DALYs = TRUE: outcomes are DALYs incurred, so fewer is better (as in IncrementalAnalysis)
genProbCE <- function(psa_output, wtp, DALYs = TRUE) {

  # psa_res is a list with elements: costs, outcomes
  costs <- as.matrix(psa_output$costs)
  outcomes <- as.matrix(psa_output$outcomes)

  # effectiveness: higher is better (DALYs incurred are negated)
  # (the Excel model, and this function until Nov 2026, added DALYs: NMB = DALYs * wtp - COSTS)
  eff <- if (DALYs) -outcomes else outcomes
  NMB <- eff*wtp - costs                        # NMB = -DALYs * wtp - COSTS
  #maxNMB <- apply(NMB, MARGIN = 1, max)         # find maximum for each row (max nmb)
  maxNMB <- max.col(NMB, ties.method = "first")  ## much faster
  #whichMaxNMB <- (NMB == maxNMB)                # check which intervention has max NMB at that WTP, it is the most cost effective at that WTP
  #return(apply(whichMaxNMB, 2, mean))           # calculate probability of cost-effectiveness
  
  # tabulate probabilities
  prob <- tabulate(maxNMB, nbins = ncol(NMB)) / nrow(NMB)
  
  # keep intervention names
  names(prob) <- colnames(NMB)
  return(prob)
}


RunCEAC <- function(psa_output, wtp, by, DALYs = TRUE) {

  wtps <- seq(0, wtp, by = by)
  
  # convert the draws to matrices once, not once per willingness-to-pay value
  psa_m <- list(costs = as.matrix(psa_output$costs), outcomes = as.matrix(psa_output$outcomes))
  
  # probability of cost effectiveness: a row per willingness-to-pay value, a column per drug
  probCE <- t(vapply(wtps, function(x) genProbCE(psa_m, wtp = x, DALYs = DALYs), numeric(ncol(psa_m$costs))))
  
  return(data.frame(WTP = wtps, probCE, check.names = FALSE))
  
}



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

## Markov trace, DALYs and costs for one treatment sequence
# primary, secondary: parameters of the first- and second-line drug, in this order:
#   c(disAE, disEff, disOther, pHosp, pWt, eps, diab, cost, qty)
# pRemit: probability of remission when treated in hospital
# p: list of model parameters plus the quantities RunModel() derives from them
#    (death probabilities, life expectancy, discount factors, costs per state, ...)
mtrace <- function(primary, secondary, pRemit, p) {
  with(p, {
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

    trans[4,6,] <- trans[5,7,] <- pHospUntrt

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
    YLL <- c(0, int.YLL[-nCycles] * expLife[-1] * pDie[-nCycles])

    # calculate vector of utilities then count the DALYs
    res.1st.line <- residualDALY + disutilWeight * primary[5] + eps.risk(pEPSPLA,primary[6]) * durEPS * disutilEPS + diab.risk(pDiabetesHAL,primary[7]) * disutilDIAB
    res.2nd.line <- residualDALY + disutilWeight * secondary[5] + eps.risk(pEPSPLA,secondary[6]) * durEPS * disutilEPS + diab.risk(pDiabetesHAL,secondary[7]) * disutilDIAB
    res.3rd.line <- residualDALY + disutilWeight * pWtFLU + eps.risk(pEPSPLA,epsFLU) * durEPS * disutilEPS + pDiabetesFLU * disutilDIAB
    res.off.1st.line <- res.off.2nd.line <- residualDALY
    ac.off.1st.line <- ac.off.2nd.line <- ac.was.on.1st.line <- ac.was.on.2nd.line <- acuteDALY + disutilWeight * pWtCPZ + eps.risk(pEPSPLA,epsCPZ) * durEPS * disutilEPS + diab.risk(pDiabetesHAL,diabCPZ) * disutilDIAB
    utils.vec <- c(res.1st.line, res.2nd.line, res.3rd.line, res.off.1st.line, res.off.2nd.line, ac.off.1st.line, ac.off.2nd.line, ac.was.on.1st.line, ac.was.on.2nd.line, 0)

    YLD <- pop.trace %*% utils.vec
    DALYS <- discFactOutcomes * (YLL + YLD)

    # calculate vector of costs and compute total costs
    c.res.1st <- op.drug.costs(primary[8], primary[9]) + eps.risk(pEPSPLA,primary[6]) * op.drug.costs(artane_cost, artane_qty) * durEPS + diab.risk(pDiabetesHAL,primary[7]) * cDiabetes + op.personnel.costs + op.lab.costs + op.capital.costs + op.carer.travel.costs + op.patient.travel.costs + op.indirect + carer.indirect
    c.res.2nd <- op.drug.costs(secondary[8], secondary[9]) + eps.risk(pEPSPLA,secondary[6]) * op.drug.costs(artane_cost, artane_qty) * durEPS + diab.risk(pDiabetesHAL,secondary[7]) * cDiabetes + op.personnel.costs + op.lab.costs + op.capital.costs + op.carer.travel.costs + op.patient.travel.costs + op.indirect + carer.indirect
    c.res.3rd <- fluphenazine_cost * fluphenazine_qty * 12 + eps.risk(pEPSPLA,epsFLU) * op.drug.costs(artane_cost, artane_qty) * durEPS + pDiabetesFLU * cDiabetes + op.personnel.costs + op.lab.costs + op.capital.costs + op.carer.travel.costs + op.patient.travel.costs + op.indirect + carer.indirect
    c.res.off.1st <- c.res.off.2nd <- 0
    c.ac.off.1st <- c.ac.off.2nd <- c.ac.was.on.1st <- c.ac.was.on.2nd <- ip.drug.costs(chlorpromazine_cost, chlorpromazine_qty, LOS) + eps.risk(pEPSPLA,epsCPZ) * ip.drug.costs(artane_cost, artane_qty, LOS) + diab.risk(pDiabetesHAL,diabCPZ) * cDiabetes + ip.personnel.costs + ip.lab.costs + ip.capital.costs + ip.carer.travel.costs + ip.patient.travel.costs + ip.indirect + carer.indirect

    cost.vec <- c(c.res.1st, c.res.2nd, c.res.3rd, c.res.off.1st, c.res.off.2nd, c.ac.off.1st, c.ac.off.2nd, c.ac.was.on.1st, c.ac.was.on.2nd, 0)

    costs <- discFactCosts * (pop.trace %*% cost.vec)

    # return objects we can work with
    #return(list(trace = pop.trace, YLL = YLL, YLD.trace = YLD.trace, YLD = YLD, 
    #            tot.DALYs = DALYS, DALYS = sum(DALYS), 
    #            cost.trace = cost.trace, cost.vector = cost.vec, tot.costs = costs, Costs = sum(costs)))

    # return objects we can work with
    # trace: proportion of the cohort in each state in each cycle
    list(DALYS = sum(DALYS), Costs = sum(costs), trace = pop.trace)
  })
}


RunModel <- function(inputs,
                     baseparms,
                     basetransitions,
                     DALYs = TRUE,
                     summary_only = FALSE) {
  
  # model parameters, used by name in the formulas below
  p <- c(as.list(baseparms), as.list(basetransitions),
         list(nStates = inputs$nStates, nCycles = inputs$nCycles))
  
  # Background mortality: the life table gives the probability of dying within each age band
  # (5-year bands in the WHO life table for Uganda), so convert it to a yearly probability
  lt <- inputs$lifeTab
  band <- c(diff(lt$ages), tail(diff(lt$ages), 1))
  q_year <- 1 - (1 - lt$p.death)^(1 / band)
  
  # add the quantities derived from the parameters
  p <- within(p, {
    cycles <- seq(0, nCycles-1, by = time)
  
    # discount factor
    discFactCosts <- 1 / (1 + discCosts) ^ cycles
    discFactOutcomes <- 1 / (1 + discOutcomes) ^ cycles
  
    ## age in each cycle
    ages <- age + cycles
  
    # match age in each cycle with a probability of death
    p.death <- q_year[match(findInterval(x=ages,vec=inputs$lifeTab$ages),inputs$lifeTab$index)]
    p.Schizo.death <- inputs$schizoTab$p.death[match(findInterval(x=ages,vec=inputs$schizoTab$ages),inputs$schizoTab$index)]
    expLife <- inputs$lifeExp$yrs[match(findInterval(x=ages,vec=inputs$lifeExp$ages),inputs$lifeExp$index)]
  
    # calculate probability of death (from other causes, schizophrenia specific -- but not both)
    pDie <- p.death + p.Schizo.death - (p.death * p.Schizo.death)
  
    # diabetes mean disutility (but put into RunModel function since they will need to be run through sensitivity analysis)
    disutilDIAB <- pNeuro * dNeuro + pFoot * dFoot + pImpotent * dImpotent + pNephro * dNephro + dDiabetes
  
    # no need for function for these below as they do not depend on the alternative chosen
    # b) Personnel costs (a complicated formula)
    ip.personnel.costs <- (((1/7)*(major_week*PSY_major*(duration_major/60)*(rate_psy/workhours_year))+(major_week*PCO_major*(duration_major/60)*(rate_pco/workhours_year))+(major_week*PNO_major*(duration_major/60)*(rate_pno/workhours_year))+(minor_week*PCO_minor*(duration_minor/60)*(rate_pco/workhours_year))+(minor_week*PNO_minor*(duration_minor/60)*(rate_pno/workhours_year)))+(PNO_drugadmin*(duration_drugadmin/60)))*time*LOS
    op.personnel.costs <- (((duration_consult/60)*(PSY*PSY_prop*(rate_psy/workhours_year))+(PCO*(1-PSY_prop)*(rate_pco/workhours_year))+(nurse*nurse_prop*(rate_pno/workhours_year)))+(pharm*(duration_drugdisp/60)*(rate_pharm/workhours_year)))*visitation_frequency*time*12
  
    # c) laboratory costs
    ip.lab.costs <- cpi_ratio(inputs$cpi, 2001)*(FBC_cost*FBC_freq+ESR_cost*ESR_freq+TPHA_cost*TPHA_freq*TPHA_proportion+VDRL_cost*VDRL_freq*(1-TPHA_proportion)+HIV_cost*HIV_freq)*time*12
    op.lab.costs <- (cpi_ratio(inputs$cpi, 2001)*(FBC_cost*FBC_freq+ESR_cost*ESR_freq))*time*12
  
    # d) capital and overheads
    ip.capital.costs <- capital_costs_per_inpatient_per_day*time*LOS
    op.capital.costs <- capital_costs_per_outpatient_visit*visitation_frequency*time*12
  
    # e) travel costs
    ip.carer.travel.costs <- averagetravel_cost*2*caregiver_prop*(time*LOS*caregiver_freq)/every
    op.carer.travel.costs <- caregiver_prop*averagetravel_cost*2*visitation_frequency*time*12
  
    ip.patient.travel.costs <- averagetravel_cost*2
    op.patient.travel.costs <- averagetravel_cost*2*visitation_frequency*time*12
  
    # f) Indirect costs
    carer.indirect <- laborforcerate*(1-unemploymentrate)*GDP_percapita*time*DALY_caretaker*(264/365)
    ip.indirect <- laborforcerate*(1-unemploymentrate)*GDP_percapita*time
    op.indirect <- laborforcerate*(1-unemploymentrate)*(GDP_percapita/264)*(visitation_frequency*time*12)
  })
  
  ## Markov traces for each treatment sequence (in the order of inputs$interventions)
  traces <- with(p, {
    mTrace.cpz <- mtrace(primary = c(disAECPZ, disEffCPZ, disOtherCPZ, pHospCPZ, pWtCPZ, epsCPZ, diabCPZ, chlorpromazine_cost, chlorpromazine_qty), 
                         secondary = c(disAERIS, disEffRIS, disOtherRIS, pHospRIS, pWtRIS, epsRIS, diabRIS, risperidone_cost, risperidone_qty), 
                         pRemit = pRemitCPZ, p = p)
  
    mTrace.hal <- mtrace(primary = c(disAEHAL, disEffHAL, disOtherHAL, pHospHAL, pWtHAL, epsHAL, diabHAL, haloperidol_cost, haloperidol_qty), 
                         secondary = c(disAERIS, disEffRIS, disOtherRIS, pHospRIS, pWtRIS, epsRIS, diabRIS, risperidone_cost, risperidone_qty), 
                         pRemit = pRemitCPZ, p = p)
  
    mTrace.ris <- mtrace(primary = c(disAERIS, disEffRIS, disOtherRIS, pHospRIS, pWtRIS, epsRIS, diabRIS, risperidone_cost, risperidone_qty),
                         secondary = c(disAEHAL, disEffHAL, disOtherHAL, pHospHAL, pWtHAL, epsHAL, diabHAL, haloperidol_cost, haloperidol_qty), 
                         pRemit = pRemitCPZ, p = p)
  
    mTrace.ola <- mtrace(primary = c(disAEOLA, disEffOLA, disOtherOLA, pHospOLA, pWtOLA, epsOLA, diabOLA, olanzapine_cost, olanzapine_qty),
                         secondary = c(disAEHAL, disEffHAL, disOtherHAL, pHospHAL, pWtHAL, epsHAL, diabHAL, haloperidol_cost, haloperidol_qty), 
                         pRemit = pRemitCPZ, p = p)
  
    mTrace.que <- mtrace(primary = c(disAEQUE, disEffQUE, disOtherQUE, pHospQUE, pWtQUE, epsQUE, diabQUE, quetiapine_cost, quetiapine_qty),
                         secondary = c(disAEHAL, disEffHAL, disOtherHAL, pHospHAL, pWtHAL, epsHAL, diabHAL, haloperidol_cost, haloperidol_qty), 
                         pRemit = pRemitCPZ, p = p)
    list(mTrace.cpz, mTrace.hal, mTrace.que, mTrace.ris, mTrace.ola)
  })
  
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
  
  # # probability of cost effectiveness calculations
  # probCE <- map_dfr(wtps, .f = genProbCE, psa_output) %>%
  #   # bind_rows() %>%
  #   bind_cols(WTP = wtps, .)
  
  probCE <- t(sapply(wtps, FUN = function(x) genProbCE(psa_output, wtp = x, DALYs = DALYs))) %>%
    as_tibble() %>%
    bind_cols(WTP = wtps, .)
  
  # use apply functions because they are easier to understand
  
  return(probCE)
  
}



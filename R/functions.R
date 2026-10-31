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
#    read from the CPI sheet of data_input.xlsx. Costs are reported in 2012 US$.
cpi_ratio <- function(from, to = 2012) {
  cpi$Value[cpi$Year == to] / cpi$Value[cpi$Year == from]
}

# 4) costs models
# a) drug costs
ip.drug.costs <- function(cost, qty) {
  return(cost * qty * 60)
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

## function to return markov traces for each intervention
## INSTRUCTIONS
# 1) first argument is a vector of parameters for first line agent in this order: c(disAE, disEff, disOther, pHosp, pWt, eps, diab, cost, qty)
# 2) second argument is a vector of parameters for second line agent in this order: c(disAE, disEff, disOther, pHosp, pWt, eps, diab, cost, qty)
# 3) third argument is the probability of hospitalization of untreated
# 4) forth argument is the probability of remission when treated in the hospital
# 5) fifth argument is the number of states (creates the two dimensional TP matrix)
# 6) sixth argument is the number of cycles (to create the third dimension that captures time varying death probabilities)
mtrace <- function(primary = numeric(), 
                   secondary = numeric(), 
                   pHospUntrt, 
                   pRemit, 
                   nStates, 
                   nCycles) {
  ## Set up 3 dimensional transition probability arrays for comparators (3rd dimension captures time varying death probabilities)
  ## Then run memory based traces
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
  
  ## get time varying death probabilities
  for (i in 1:nCycles) {
    # insert death probabilities
    trans[,nStates,i] <- pDie[i]
    # insert probability of remaining in the same state in each cycle (complement of rowsums)
    diag(trans[,,i]) <- 1-rowSums(trans[,,i])
  }
  # set absorbing state (death) probability to 1
  trans[nStates,nStates,] <- 1
  
  # NOW RUN THE MODEL TO GENERATE TRACES 
  #Initialize population trace matrices [cycle,state,sojourn time]
  trace <- array(0,c(nCycles,nStates,nCycles+1))
  trace[1,,1] <- c(1,0,0,0,0,0,0,0,0,0)
  
  # index the cycles by c, and for each cycle calculate the proportions transiting 
  for (c in 2:nCycles) {
    for (cs in 1:(c-1)) {
      ## get cycle specific transition matrix
      cond.trans <- trans[,,cs]
      # how many people transition from cycle c-1 to c? a matrix multiplication
      new.state <- trace[c-1,,cs] %*% cond.trans
      # update the trace array (keeping track of time in state)
      trace[c,,cs+1] <- new.state + trace[c,,cs]
    }
  }
  
  # add up time spend in each state
  pop.trace <- apply(trace,c(1,2),sum)
  
  # calculate YLL
  int.YLL <- rowSums(pop.trace[,1:9])
  YLL <- numeric()
  YLL[1] <- 0
  for (i in 2:nCycles) {
    YLL[i] <- int.YLL[i-1] * expLife[i] * pDie[i-1]
  }
  
  # calculate vector of utilities then count the DALYs
  res.1st.line <- residualDALY + disutilWeight * primary[5] + eps.risk(pEPSPLA,primary[6]) * durEPS * disutilEPS + diab.risk(pDiabetesHAL,primary[7]) * disutilDIAB
  res.2nd.line <- residualDALY + disutilWeight * secondary[5] + eps.risk(pEPSPLA,secondary[6]) * durEPS * disutilEPS + diab.risk(pDiabetesHAL,secondary[7]) * disutilDIAB
  res.3rd.line <- residualDALY + disutilWeight * pWtFLU + eps.risk(pEPSPLA,epsFLU) * durEPS * disutilEPS + pDiabetesFLU * disutilDIAB
  res.off.1st.line <- res.off.2nd.line <- residualDALY
  ac.off.1st.line <- ac.off.2nd.line <- ac.was.on.1st.line <- ac.was.on.2nd.line <- acuteDALY + disutilWeight * pWtCPZ + eps.risk(pEPSPLA,epsCPZ) * durEPS * disutilEPS + diab.risk(pDiabetesHAL,diabCPZ) * disutilDIAB
  utils.vec <- c(res.1st.line, res.2nd.line, res.3rd.line, res.off.1st.line, res.off.2nd.line, ac.off.1st.line, ac.off.2nd.line, ac.was.on.1st.line, ac.was.on.2nd.line, 0)
  
  YLD <- pop.trace %*% utils.vec
  DALYS <- discFactOutcomes * (YLL + YLD)
  
  YLD.trace <- cbind(pop.trace[,1] * res.1st.line, pop.trace[,2] * res.2nd.line, pop.trace[,3] * res.3rd.line, pop.trace[,4] * res.off.1st.line, pop.trace[,5] * res.off.2nd.line, pop.trace[,6] * ac.off.1st.line, pop.trace[,7] * ac.off.1st.line, pop.trace[,8] * ac.off.1st.line, pop.trace[,9] * ac.off.1st.line)
  
  # calculate vector of costs and compute total costs
  c.res.1st <- op.drug.costs(primary[8], primary[9]) + eps.risk(pEPSPLA,primary[6]) * ip.drug.costs(artane_cost, artane_qty) + diab.risk(pDiabetesHAL,primary[7]) * cDiabetes + op.personnel.costs + op.lab.costs + op.capital.costs + op.carer.travel.costs + op.patient.travel.costs + op.indirect + carer.indirect
  c.res.2nd <- op.drug.costs(secondary[8], secondary[9]) + eps.risk(pEPSPLA,secondary[6]) * ip.drug.costs(artane_cost, artane_qty) + diab.risk(pDiabetesHAL,secondary[7]) * cDiabetes + op.personnel.costs + op.lab.costs + op.capital.costs + op.carer.travel.costs + op.patient.travel.costs + op.indirect + carer.indirect
  c.res.3rd <- fluphenazine_cost * fluphenazine_qty * 12 + eps.risk(pEPSPLA,epsFLU) * ip.drug.costs(artane_cost, artane_qty) + pDiabetesFLU * cDiabetes + op.personnel.costs + op.lab.costs + op.capital.costs + op.carer.travel.costs + op.patient.travel.costs + op.indirect + carer.indirect
  c.res.off.1st <- c.res.off.2nd <- 0
  # Artane price index 145.14 is not in the CPI sheet (2009 = 144.24); kept as in the 2015 Excel model until the source year is confirmed
  c.ac.off.1st <- c.ac.off.2nd <- c.ac.was.on.1st <- c.ac.was.on.2nd <- ip.drug.costs(chlorpromazine_cost, chlorpromazine_qty)*cpi_ratio(2010) + eps.risk(pEPSPLA,epsCPZ) * ip.drug.costs(artane_cost, artane_qty)*(196.46/145.14) + diab.risk(pDiabetesHAL,diabCPZ) * cDiabetes + ip.personnel.costs + ip.lab.costs + ip.capital.costs + ip.carer.travel.costs + ip.patient.travel.costs + ip.indirect + carer.indirect
  
  cost.trace <- cbind(pop.trace[,1] * c.res.1st, pop.trace[,2] * c.res.2nd, pop.trace[,3] * c.res.3rd, pop.trace[,4] * c.res.off.1st, pop.trace[,5] * c.res.off.2nd, pop.trace[,6] * c.ac.off.1st, pop.trace[,7] * c.ac.off.2nd, pop.trace[,8] * c.ac.was.on.1st, pop.trace[,9] * c.ac.was.on.2nd)
  
  cost.vec <- c(c.res.1st, c.res.2nd, c.res.3rd, c.res.off.1st, c.res.off.2nd, c.ac.off.1st, c.ac.off.2nd, c.ac.was.on.1st, c.ac.was.on.2nd, 0)
  
  costs <- discFactCosts * (pop.trace %*% cost.vec)
  
  # return objects we can work with
  #return(list(trace = pop.trace, YLL = YLL, YLD.trace = YLD.trace, YLD = YLD, 
  #            tot.DALYs = DALYS, DALYS = sum(DALYS), 
  #            cost.trace = cost.trace, cost.vector = cost.vec, tot.costs = costs, Costs = sum(costs)))
  
  # return objects we can work with
  return(list(DALYS = sum(DALYS), Costs = sum(costs)))
  
}


RunModel <- function(interventions,
                     baseparms,
                     basetransitions,
                     DALYs = TRUE,
                     summary_only = FALSE) {
  
  # need variable names from each element of these named vectors (improve this later)
  parms <- c(baseparms, basetransitions)
  list2env(as.list(parms), envir = environment())  # creates named parameters in *this* environment
  #beta1 <- as.list(parms)
  #attach(beta1)
  
  
  # ---------------------------------------------------------------------------------
  # makes temp copy of mtrace function so variables created here are accessible to it
  # without needing to place entire mtrace function in the RunModel function
  # ---------------------------------------------------------------------------------
  
  environment(mtrace) <- environment() 

  cycles <- seq(0, nCycles-1, by = time)
  
  # discount factor
  discFactCosts <- 1 / (1 + discCosts) ^ cycles
  discFactOutcomes <- 1 / (1 + discOutcomes) ^ cycles
  
  ## age in each cycle
  ages <- age + cycles
  
  # match age in each cycle with a probability of death
  p.death <- lifeTab$p.death[match(findInterval(x=ages,vec=lifeTab$ages),lifeTab$index)]
  p.Schizo.death <- schizoTab$p.death[match(findInterval(x=ages,vec=schizoTab$ages),schizoTab$index)]
  expLife <- lifeExp$yrs[match(findInterval(x=ages,vec=lifeExp$ages),lifeExp$index)]
  
  # calculate probability of death (from other causes, schizophrenia specific -- but not both)
  pDie <- p.death + p.Schizo.death - (p.death * p.Schizo.death)
  
  # diabetes mean disutility (but put into RunModel function since they will need to be run through sensitivity analysis)
  disutilDIAB <- pNeuro * dNeuro + pFoot * dFoot + pImpotent * dImpotent + pNephro * dNephro + dDiabetes
  
  # no need for function for these below as they do not depend on the alternative chosen
  # b) Personnel costs (a complicated formula)
  ip.personnel.costs <- (((1/7)*(major_week*PSY_major*(duration_major/60)*(rate_psy/workhours_year))+(major_week*PCO_major*(duration_major/60)*(rate_pco/workhours_year))+(major_week*PNO_major*(duration_major/60)*(rate_pno/workhours_year))+(minor_week*PCO_minor*(duration_minor/60)*(rate_pco/workhours_year))+(minor_week*PNO_minor*(duration_minor/60)*(rate_pno/workhours_year)))+(PNO_drugadmin*(duration_drugadmin/60)))*time*LOS
  op.personnel.costs <- (((duration_consult/60)*(PSY*PSY_prop*(rate_psy/workhours_year))+(PCO*(1-PSY_prop)*(rate_pco/workhours_year))+(nurse*nurse_prop*(rate_pno/workhours_year)))+(pharm*(duration_drugdisp/60)*(rate_pharm/workhours_year)))*visitation_frequency*time*12
  
  # c) laboratory costs
  ip.lab.costs <- cpi_ratio(2001)*(FBC_cost*FBC_freq+ESR_cost*ESR_freq+TPHA_cost*TPHA_freq*TPHA_proportion+VDRL_cost*VDRL_freq*(1-TPHA_proportion)+HIV_cost*HIV_freq)*time*12
  op.lab.costs <- (cpi_ratio(2001)*(FBC_cost*FBC_freq+ESR_cost*ESR_freq))*time*12
  
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
  
  
  ## apply mtrace function to each of the drugs to return markov traces for each intervention
  
  mTrace.cpz <- mtrace(primary = c(disAECPZ, disEffCPZ, disOtherCPZ, pHospCPZ, pWtCPZ, epsCPZ, diabCPZ, chlorpromazine_cost, chlorpromazine_qty), 
                       secondary = c(disAERIS, disEffRIS, disOtherRIS, pHospRIS, pWtRIS, epsRIS, diabRIS, risperidone_cost, risperidone_qty), 
                       pHospUntrt = pHospUntrt, 
                       pRemit = pRemitCPZ,
                       nStates = nStates, 
                       nCycles = nCycles)
  
  mTrace.hal <- mtrace(primary = c(disAEHAL, disEffHAL, disOtherHAL, pHospHAL, pWtHAL, epsHAL, diabHAL, haloperidol_cost, haloperidol_qty), 
                       secondary = c(disAERIS, disEffRIS, disOtherRIS, pHospRIS, pWtRIS, epsRIS, diabRIS, risperidone_cost, risperidone_qty), 
                       pHospUntrt = pHospUntrt, 
                       pRemit = pRemitCPZ,
                       nStates = nStates, 
                       nCycles = nCycles)
  
  mTrace.ris <- mtrace(primary = c(disAERIS, disEffRIS, disOtherRIS, pHospRIS, pWtRIS, epsRIS, diabRIS, risperidone_cost, risperidone_qty),
                       secondary = c(disAEHAL, disEffHAL, disOtherHAL, pHospHAL, pWtHAL, epsHAL, diabHAL, haloperidol_cost, haloperidol_qty), 
                       pHospUntrt = pHospUntrt, 
                       pRemit = pRemitCPZ,
                       nStates = nStates, 
                       nCycles = nCycles)
  
  mTrace.ola <- mtrace(primary = c(disAEOLA, disEffOLA, disOtherOLA, pHospOLA, pWtOLA, epsOLA, diabOLA, olanzapine_cost, olanzapine_qty),
                       secondary = c(disAEHAL, disEffHAL, disOtherHAL, pHospHAL, pWtHAL, epsHAL, diabHAL, haloperidol_cost, haloperidol_qty), 
                       pHospUntrt = pHospUntrt, 
                       pRemit = pRemitCPZ,
                       nStates = nStates, 
                       nCycles = nCycles)
  
  mTrace.que <- mtrace(primary = c(disAEQUE, disEffQUE, disOtherQUE, pHospQUE, pWtQUE, epsQUE, diabQUE, quetiapine_cost, quetiapine_qty),
                       secondary = c(disAEHAL, disEffHAL, disOtherHAL, pHospHAL, pWtHAL, epsHAL, diabHAL, haloperidol_cost, haloperidol_qty), 
                       pHospUntrt = pHospUntrt, 
                       pRemit = pRemitCPZ,
                       nStates = nStates, 
                       nCycles = nCycles)
  
  # don't forget to detach
  #detach(beta1)
  
  #interventions <- c("Chlorpromazine", "Haloperidol", "Quetiapine", "Risperidone", "Olanzapine")
  outcomes <- c(mTrace.cpz$DALYS, mTrace.hal$DALYS, mTrace.que$DALYS, mTrace.ris$DALYS, mTrace.ola$DALYS)
  costs <- c(mTrace.cpz$Costs, mTrace.hal$Costs, mTrace.que$Costs, mTrace.ris$Costs, mTrace.ola$Costs)
  
  output_tab <- tibble(intervention = interventions,
                       costs = costs,
                       outcomes = outcomes)
  
  # establish efficiency frontier (CEEF)
  # pass outputs into efficiency frontier function
  # ceef_analysis <- CE_Frontier(intervention = interventions, costs = costs, outcomes = outcomes, DALYs = TRUE)
  
  # return things in output I need for other functions
  if(summary_only) return(list(costs=costs, outcomes=outcomes)) else 
    return(list(ce_result = output_tab, 
              baseparms = baseparms, 
              basetransitions = basetransitions, 
              DALYs = DALYs))
  
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


## loop through base-case, high and low values and calculate model
# function requires user to have run the incremental analysis, if not, ---> stop, ask user to run basecase first
# run incremental within the function and extract the things we need
owsa <- function(model, 
                 low_base, low_transitions, 
                 high_base, high_transitions, 
                 outcome, 
                 intervention,
                 comparator) {

  # Result of interest: `intervention` versus `comparator` --
  # incremental cost ("cost"), DALYs averted / QALYs gained ("outcome") or ICER ("ICER")
  pairwise <- function(m) {
    r  <- m$ce_result
    dc <- r$costs[r$intervention == intervention] - r$costs[r$intervention == comparator]
    do <- if (m$DALYs) r$outcomes[r$intervention == comparator] - r$outcomes[r$intervention == intervention]
          else r$outcomes[r$intervention == intervention] - r$outcomes[r$intervention == comparator]
    switch(outcome, ICER = dc / do, cost = dc, outcome = do)
  }
  
  # generate base case from model object
  interventions <- model$ce_result %>%
    pull(intervention)
  base_case <- RunModel(interventions = interventions,
                        baseparms = model$baseparms, 
                        basetransitions = model$basetransitions,
                        DALYs = model$DALYs)
  
  base.case <- pairwise(base_case)
  
  # set up matrix to collect the results
  owsamat <- matrix(NA,
                    nrow = length(low_base) + length(low_transitions),
                    ncol=2,
                    dimnames=list(c(names(low_base), names(low_transitions)),
                                  c("low", "high")))
  
  lowinputs <- low_base
  highinputs <- high_base
  owsavars <- c(paste0(names(low_base), " (", sprintf(lowinputs, fmt = "%#.2f"),", ", sprintf(highinputs, fmt = "%#.2f"),")"),
                paste0("TP set - ", names(low_transitions)))
  
  i=0
  # set all live values to the base case
  live <- model$baseparms
  for (n in names(lowinputs)) {
    i=i+1
    # counter
    message("Running DSA on ", n, ", variable ", i, " of ", length(names(lowinputs)))
    # replace the current parameter with its low value and run the model
    low.input <- lowinputs[n]
    live[n] <- low.input
    dsa.low <- RunModel(interventions = interventions,
                        baseparms = live, 
                        basetransitions = model$basetransitions)
    low_val <- pairwise(dsa.low)
    
    # replace the current parameter with its high value and run the model
    high.input <- highinputs[n]
    live[n] <- high.input
    dsa.high <- RunModel(interventions = interventions,
                         baseparms = live, 
                         basetransitions = model$basetransitions)
    high_val <- pairwise(dsa.high)
    
    owsamat[i,] <- c(low_val, high_val)

    # reset to the original value
    live[n] <- model$baseparms[n]
  }
  
  # do the same for the transition probabilities, but put all in at once
  # set all live values to the base case, the 4th element of list returned by the RunModel function
  live <- model$basetransitions
  
  for (i in 1:length(low_transitions)) {

    message("Running DSA on TP set for ", names(low_transitions)[i])
    
    # replace the current parameter-set with its low value and run the model
    low.input <- low_transitions[[i]]
    live[names(low_transitions[[i]])] <- low.input
    
    dsa.low <- RunModel(interventions = interventions,
                        baseparms = model$baseparms, 
                        basetransitions = live)
    low_val <- pairwise(dsa.low)
    
    # replace the current parameter-set with its high value and run the model
    high.input <- high_transitions[[i]]
    live[names(high_transitions[[i]])] <- high.input
    dsa.high <- RunModel(interventions = interventions,
                         baseparms = model$baseparms, 
                         basetransitions = live)
    high_val <- pairwise(dsa.high)

    owsamat[i+length(low_base), ] <- c(low_val, high_val)
    
    # reset to the original value in basetransitions
    live[names(high_transitions[[i]])] <- model$basetransitions[names(high_transitions[[i]])]
    
  }
  
  # post processing of data to return full results
  owsa_tab <- as_tibble(owsamat) %>%
    add_column(variable = owsavars, .before = "low") %>%
    mutate(range = abs(high - low), high_val = high - base.case, low_val = low - base.case) %>%
    arrange(range)
  
  return(list(owsa_data = owsa_tab, owsa_outcome = outcome, owsa_basecase = base.case, DALYs = model$DALYs))
  
}

# function takes output from DSA as a list, and plots the results (user inputs max vars)
TornadoPlot <- function(owsa_list, max_vars) {
  
  # extract the things we need
  owsa_tab <- owsa_list$owsa_data
  outcome <- owsa_list$owsa_outcome
  basecase <- owsa_list$owsa_basecase
  DALYs <- owsa_list$DALYs
  on_outcome <- ifelse(DALYs == TRUE, "DALYs averted", "QALYs gained")
  
  # set color
  col <- brewer.pal(3,"Dark2")
  
  # filter for those with perc impact on outcome
  owsa_torn <- owsa_tab %>%
    slice_max(range, n = max_vars) %>%
    arrange(range)
  
  axislabel <- switch(outcome,
                      ICER=paste0("Cost per ", on_outcome),
                      cost="Incremental Cost ($)",
                      outcome=on_outcome)
  
  plot_ly(data = owsa_torn) %>%
    add_trace(x = ~low_val, y = ~variable, 
              base = basecase, 
              type = 'bar', orientation = 'h', name = "lower input",
              marker = list(color = "rgba(166,206,227,0.7)",
                            line = list(color = "rgb(166,206,227)",
                                        width = 2))) %>%
    add_trace(x = ~high_val, y = ~variable, 
              base = basecase, 
              type = 'bar', orientation = 'h', name = "higher input",
              marker = list(color = "rgba(31,120,180,0.7)",
                            line = list(color = "rgb(31,120,180)",
                                        width = 2))) %>%
    layout(barmode='overlay') %>% 
    layout(yaxis = list(categoryorder = "array",
                        categoryarray = ~variable,
                        title = ""),
           xaxis = list(zeroline = FALSE,
                        showline = FALSE,
                        showgrid = FALSE,
                        title = axislabel),
           font = list(size=12, 
                       family='Fira Sans'))
  
}


## ----------------------------
## PSA utilities
## ----------------------------

# functions to create intervals for the dirichlet distributions
dir_interval <- function(prior) {
  prior <- unlist(prior)
  posterior <- prior + rep(1, length=length(prior))
  lower <- qgamma(.025,shape=posterior,rate=1)/sum(qgamma(.025,shape=posterior,rate=1))
  upper <- qgamma(.975,shape=posterior,rate=1)/sum(qgamma(.975,shape=posterior,rate=1))
  return(list(lower,upper))
}

# function to generate and name random dirichlet variates
vdirichlet <- function(n, alpha) {
  variates <- MCMCpack::rdirichlet(n, alpha)
  colnames(variates) <- names(alpha)
  return(variates)
}

## process draws inside function for improved readability -- this function is only called once per analysis run 
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
  
  # ... then exactly the same code for Beta/Normal/Lognormal/Dirichlet draws ...
  ## use tidyverse styling for readability
  # --- Beta draws
  beta_dt <- map_dfc(betaparms, ~{
    tibble(!!.x := rbeta(nsims, alphas[.x], betas[.x]))
  })
  
  # --- Normal draws
  normal_dt <- map_dfc(normparms, ~{
    tibble(!!.x := rnorm(nsims, mean = base_vals[.x], sd = sems[.x]))
  })
  
  # --- Lognormal draws
  logn_dt <- map_dfc(lognormparms, ~{
    mu    <- log(base_vals[.x])
    sigma <- (log(high_vals[.x]) - log(low_vals[.x])) / 3.92
    tibble(!!.x := rlnorm(nsims, meanlog = mu, sdlog = sigma))
  })
  
  # --- Dirichlet draws
  dir_dt <- {
    dir_rows <- transitions %>%
      filter(dist == "dirichlet") %>%
      pull(alpha, variable)
    
    randdirchinputs <- dir_rows %>%
      split(., cut(seq_along(.),
                   length(.)/length(interventions),
                   labels = interventions))
    
    map(randdirchinputs, function(alpha) {
      k <- length(alpha)
      G <- matrix(rgamma(nsims * k, shape = rep(alpha, each = nsims)),
                  nrow = nsims)
      D <- G / rowSums(G)
      colnames(D) <- names(alpha)
      as_tibble(D)
    }) %>% bind_cols()
  }
  
  # --- Fixed parameters
  fixed_dt <- if (length(fixedparms) > 0) {
    tibble(!!!map(set_names(fixedparms), ~ rep(base_vals[.x], nsims)))
  } else tibble()
  
  # --- Combine all blocks
  return(bind_cols(beta_dt, normal_dt, logn_dt, dir_dt, fixed_dt))
  
}

# this is a wrapper function that calls the RunModel function on psalist
# use it to extract the PSA trace of costs and outcomes
genPSA <- function(model, psalist) {
  
  # get things we need
  basenames <- names(model$baseparms)
  transitnames <- names(model$basetransitions)
  
  interventions <- model$ce_result[["intervention"]]
  
  res <- RunModel(interventions = interventions,
                  baseparms = psalist[basenames],
                  basetransitions = psalist[transitnames])
  
  costs <- setNames(res$ce_result[["costs"]], res$ce_result[["intervention"]])
  outcomes <- setNames(res$ce_result[["outcomes"]], res$ce_result[["intervention"]])
  
  #rm(res); gc()
  
  return(list(costs = costs, outcomes = outcomes))
  
}

# psa_res <- lapply(split(psa_input, 1:nrow(psa_input)), FUN = RunPSA, base_case = model) %>%
#   pmap(., bind_rows)

RunPSA <- function(model, nsims, parameters, transitions, interventions) {
  # --- 1. Precompute all draws once, outside of genPSA
  psalist <- draw_psa_df(
    nsims        = nsims,
    parameters   = parameters,
    transitions  = transitions,
    interventions = interventions
  )
  
  # --- 2. Split rows into list of 1-row tibbles
  rowlist <- split(psalist, seq_len(nrow(psalist)))
  psaTraces <- pblapply(rowlist, function(x) genPSA(model, x))
  
  # keep psalist intact and iterate by row index
  # psaTraces <- pblapply(seq_len(nrow(psalist)), function(i) genPSA(model, psalist[i, ]))
  
  # --- 4. Extract and bind results
  costs    <- purrr::map_dfr(psaTraces, "costs")
  outcomes <- purrr::map_dfr(psaTraces, "outcomes")
  
  return(list(costs = costs, outcomes = outcomes))
}

# function takes the maxwtp, a vector or dataframe of costs and outcomes for all interventions
# and returns probability that a given intervention is cost effective 
genProbCE <- function(psa_output, wtp) {
  
  # psa_res is a list with elements: costs, outcomes
  costs <- as.matrix(psa_output$costs)
  outcomes <- as.matrix(psa_output$outcomes)

  NMB <- outcomes*wtp - costs                   # NMB = DALYs * wtp - COSTS
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


RunCEAC <- function(psa_output, wtp, by) {
  
  wtps <- seq(0, wtp, by = by)
  
  # # probability of cost effectiveness calculations
  # probCE <- map_dfr(wtps, .f = genProbCE, psa_output) %>%
  #   # bind_rows() %>%
  #   bind_cols(WTP = wtps, .)
  
  probCE <- t(sapply(wtps, FUN = function(x) genProbCE(psa_output, wtp = x))) %>% 
    as_tibble() %>%
    bind_cols(WTP = wtps, .)
  
  # use apply functions because they are easier to understand
  
  return(probCE)
  
}



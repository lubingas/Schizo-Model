# switch off file writing if in use
if(sink.number()>0) sink()

rm(list=ls())

# define function for detaching variables
detachAllData <- function()  {
  pos.to.detach <- (1:length(search()))[substring(search(), 
                                                  first = 1, last = 8) != "package:" & search() != ".GlobalEnv" & 
                                          search() != "Autoloads" & search() != "CheckExEnv" & search() != "tools:rstudio" & search() != "TempEnv"]
  for (i in 1:length(pos.to.detach)) {
    if (length(pos.to.detach) > 0) {
      detach(pos = pos.to.detach[1])
      pos.to.detach <- (1:length(search()))[substring(search(), 
                                                      first = 1, last = 8) != "package:" & search() != 
                                              ".GlobalEnv" & search() != "Autoloads" & search() != 
                                              "CheckExEnv" & search() != "tools:rstudio" & 
                                              search() != "TempEnv"]
    }
  }
}

# run detach function
detachAllData()



library(tidyverse)
library(readxl)
library(plotly)
library(RColorBrewer)


## some subroutines for calculating various stuff
# 1) extrapyramidal side effects
eps.risk <- function(base.p, drug.OR){
  return(drug.OR * base.p/(1 - base.p + (drug.OR * base.p)))
}

# 2) Diabetes risk
diab.risk <- function(base.p, drug.RR) {
  return(base.p * drug.RR)
}


# 3) costs models
# a) drug costs
ip.drug.costs <- function(cost, qty) {
  return(cost * qty * 60)
}

op.drug.costs <- function(cost, qty) {
  return(cost * qty * 365)
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
  c.ac.off.1st <- c.ac.off.2nd <- c.ac.was.on.1st <- c.ac.was.on.2nd <- ip.drug.costs(chlorpromazine_cost, chlorpromazine_qty)*196.46/149.97 + eps.risk(pEPSPLA,epsCPZ) * ip.drug.costs(artane_cost, artane_qty)*(196.46/145.14) + diab.risk(pDiabetesHAL,diabCPZ) * cDiabetes + ip.personnel.costs + ip.lab.costs + ip.capital.costs + ip.carer.travel.costs + ip.patient.travel.costs + ip.indirect + carer.indirect
  
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
                     DALYs = TRUE) {# inputs vector houses all our parameter inputs
  
  environment(mtrace) <- environment() # makes temp copy of mtrace function so variables created here are accessible to it without needing to plance entire mtrace function in the RunModel function

  parms <- c(baseparms, basetransitions)
  beta1 <- as.list(parms)
  attach(beta1)
  
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
  
  # duabetes mean disutility (but put into RunModel function since they will need to be run through sensitivity analysis)
  disutilDIAB <- pNeuro * dNeuro + pFoot * dFoot + pImpotent * dImpotent + pNephro * dNephro + dDiabetes
  
  # no need for function for these below as they do not depend on the alternative chosen
  # b) Personnel costs (a complicated formula)
  ip.personnel.costs <- (((1/7)*(major_week*PSY_major*(duration_major/60)*(rate_psy/workhours_year))+(major_week*PCO_major*(duration_major/60)*(rate_pco/workhours_year))+(major_week*PNO_major*(duration_major/60)*(rate_pno/workhours_year))+(minor_week*PCO_minor*(duration_minor/60)*(rate_pco/workhours_year))+(minor_week*PNO_minor*(duration_minor/60)*(rate_pno/workhours_year)))+(PNO_drugadmin*(duration_drugadmin/60)))*time*LOS
  op.personnel.costs <- (((duration_consult/60)*(PSY*PSY_prop*(rate_psy/workhours_year))+(PCO*(1-PSY_prop)*(rate_pco/workhours_year))+(nurse*nurse_prop*(rate_pno/workhours_year)))+(pharm*(duration_drugdisp/60)*(rate_pharm/workhours_year)))*visitation_frequency*time*12
  
  # c) laboratory costs
  ip.lab.costs <- (196.64/82.04)*(FBC_cost*FBC_freq+ESR_cost*ESR_freq+TPHA_cost*TPHA_freq*TPHA_proportion+VDRL_cost*VDRL_freq*(1-TPHA_proportion)*HIV_freq)*time*12
  op.lab.costs <- ((196.64/82.04)*(FBC_cost*FBC_freq+ESR_cost*ESR_freq))*time*12
  
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
  detach(beta1)
  
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
  
  # gather things we need from model object
  res_tab <- model$ce_result
  costs <- res_tab %>%
    pull(costs, intervention)
  outcomes <- res_tab %>%
    pull(outcomes, intervention)
  intervention <- res_tab %>%
    pull(intervention)
  DALYs <- model$DALYs
  
  # initialize a vector where we will keep interventions on the frontier
  frontier <- c()
  
  # make a table we can use
  # 1. order interventions by DALYs in descending order (intervention with highest DALYs first)
  # 2. calculate ICER relative to intervention with highest DALYs
  # 3. determine whether ICER or dominant with highest DALY intervention as reference
  output <- res_tab %>%
    arrange(if (DALYs) desc(outcomes) else (outcomes)) %>%# arrange by outcomes, first row is first intervention on frontier
    mutate(delta_outcomes = if (DALYs) (first(outcomes) - outcomes) else (outcomes - first(outcomes)), # calculate incrementals
           delta_costs = costs - first(costs), 
           iCER = (delta_costs/delta_outcomes), # calculate iCERs relative to intervention with highests DALYs
           dom = if_else(is.nan(iCER), "Ref", if_else(iCER > 0, "ICER", "Dominant"))) # check dominance
  
  # temporarily store output so I can access at the end to select original estimates of interventions on the CEEF
  temp_output <- output %>%
    select(-c(iCER, dom))
  
  # after arrange, first row is the first element of the CEEF
  # 0.intervention 1 is the origin of the CEEF (status quo)
  
  i <- 1
  
  frontier[i] <- output$intervention[1]
  
  # loop through the possibilities
  while (nrow(output) > 2) { # stop when there are only two interventions
    # CEEF logic
    #    check on second intervention: if dominant over intervention 1, then it is next point on the CEEF,
    #    else if intervention 2 has a positive ICER and intervention 3 is dominant over intervention 1, 
    #    then drop intervention 2 (it is extended dominated)
    
    # check on intervention 2: if dominant, then it is next point on the CEEF,
    if(output$dom[2] == "Dominant") {
      
      # add this intervention to the frontier
      frontier[i+1] <- output$intervention[2]
      
      # and remove intervention 1 from the output dataframe, such that this one becomes the new starting point for the CEEF
      output <- output[-1, ]
      
    } else if (output$dom[2] =="ICER") { # if in 2nd row we have +ve iCER
      
      ## nested if to check on 3rd row
      if (output$iCER[3] > output$iCER[2]) { # if iCER for intervention in third row is greater than that in row 2
        # output$intervention[2] becomes the new origin on the CEEF and it is added to the frontier
        frontier[i+1] <- output$intervention[2]
        # and remove intervention 1 from the output dataframe
        output <- output[-1, ]
        
      } else if (output$iCER[3] < output$iCER[2]) {
        # output$intervention[2] is extended dominated, so remove it from the output dataframe (or collect in separate list)
        # output$intervention[3] becomes the new starting point on the frontier (i.e. delete both interventions 1 and 2)
        frontier[i+1] <- output$intervention[3]
        output <- output[-c(1,2), ]
        
      } else {
        stop("error")
      }
    }
    
    # after checks, identifying new point on CEEF, arrange and calculate everything again
    output <- output %>%
      arrange(if (DALYs) desc(outcomes) else (outcomes)) %>%# arrange by outcomes, first row is first intervention on frontier
      mutate(delta_outcomes = if (DALYs) (first(outcomes) - outcomes) else (outcomes - first(outcomes)), # calculate incrementals
             delta_costs = costs - first(costs), 
             iCER = (delta_costs/delta_outcomes), # calculate iCERs relative to intervention with highest
             dom = if_else(is.nan(iCER), "Ref", if_else(iCER > 0, "ICER", "Dominant"))) # check dominance
    
    i <- i + 1
    
  }
  
  # send final element to frontier vector
  frontier[i+1] <- output$intervention[2]
  
  # which interventions are not on the frontier
  extended_dom <- intervention[!(intervention %in% frontier)]
  
  on_frontier <- temp_output %>%
    filter(intervention %in% frontier)
  
  # this table will summarize the incremental outcomes we are interested in along the CEEF
  res_table <- temp_output %>%
    filter(intervention %in% frontier) %>%
    mutate(delta_costs = costs - lag(costs),
           delta_outcomes = if (DALYs) -1*(outcomes - lag(outcomes)) else (outcomes - lag(outcomes)),
           iCER = delta_costs/delta_outcomes) %>% 
    replace(is.na(.), 0) %>% 
    mutate(explanation = case_when(iCER == 0 ~ "Ref",
                                   delta_costs < 0 & delta_outcomes < 0 ~ "Dominant",
                                   delta_costs < 0 & delta_outcomes > 0 ~ paste0("Dominant over ", lag(intervention)),
                                   delta_costs > 0 & delta_outcomes > 0 ~ paste0("iCER versus ", lag(intervention)))) %>% 
    select(intervention, costs, outcomes, delta_costs, delta_outcomes, iCER, explanation)
  
  return(list(ce_points = temp_output, 
              ce_lines = on_frontier, 
              incrementals = res_table, 
              extended_dominated = extended_dom, 
              outcome_to_plot = if(DALYs) "DALYs averted" else "QALYs gained"))
  
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
# function requires user to have run the incremental analysis, if not, run incremental within the function and extrace the things we need
owsa <- function(model, 
                 low_base, low_transitions, 
                 high_base, high_transitions, 
                 outcome, comparator) {
  
  # generate base case from model object
  interventions <- model$ce_result %>%
    pull(intervention)
  base_case <- RunModel(interventions = interventions,
                       baseparms = model$baseparms, 
                       basetransitions = model$basetransitions, 
                       DALYs = model$DALYs)
  
  frontier <- IncrementalAnalysis(model = base_case)
  
  # on what outcome, should DSA be run?
  #on_outcome <- model[[1]][[5]]

  base.case <- frontier$incrementals %>%
    filter(intervention %in% comparator) %>%
    pull(switch(outcome,
                ICER=iCER,
                cost=delta_costs,
                outcome=delta_outcomes))
  
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
    frontier.low <- IncrementalAnalysis(model = dsa.low)
    
    # extract value for outcome and comparator of interest
    low_val <- frontier.low$incrementals %>%
      filter(intervention %in% comparator) %>%
      pull(switch(outcome,
                  ICER=iCER,
                  cost=delta_costs,
                  outcome=delta_outcomes))
    
    # replace the current parameter with its high value and run the model
    high.input <- highinputs[n]
    live[n] <- high.input
    dsa.high <- RunModel(interventions = interventions,
                         baseparms = live, 
                         basetransitions = model$basetransitions)
    frontier.high <- IncrementalAnalysis(model = dsa.high)
    
    high_val <- frontier.high$incrementals %>%
      filter(intervention %in% comparator) %>%
      pull(switch(outcome,
                  ICER=iCER,
                  cost=delta_costs,
                  outcome=delta_outcomes))
    
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
    frontier.low <- IncrementalAnalysis(model = dsa.low)
    
    # extract value for outcome and comparator of interest
    low_val <- frontier.low$incrementals %>%
      filter(intervention %in% comparator) %>%
      pull(switch(outcome,
                  ICER=iCER,
                  cost=delta_costs,
                  outcome=delta_outcomes))
    
    # replace the current parameter-set with its high value and run the model
    high.input <- high_transitions[[i]]
    live[names(high_transitions[[i]])] <- high.input
    dsa.high <- RunModel(interventions = interventions,
                         baseparms = model$baseparms, 
                         basetransitions = live)
    frontier.high <- IncrementalAnalysis(model = dsa.high)
    
    # get the result we want and place it in our results matrix
    high_val <- frontier.high$incrementals %>%
      filter(intervention %in% comparator) %>%
      pull(switch(outcome,
                  ICER=iCER,
                  cost=delta_costs,
                  outcome=delta_outcomes))

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
  owsa_torn <- dsa_res$owsa_data %>%
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


# this function calls the RunModel function on psalist
# use it to extract the PSA trace of costs and outcomes
genPSA <- function(model, psalist) {
  
  # get things we need
  basenames <- names(model$baseparms)
  transitnames <- names(model$basetransitions)
  
  costs <- c()
  outcomes <- c()
  
  interventions <- model$ce_result %>%
    pull(intervention)
  
  res <- RunModel(interventions = interventions,
                  baseparms = psalist[basenames],
                  basetransitions = psalist[transitnames])
  
  costs <- res$ce_result %>%
    pull(costs, intervention)
  
  outcomes <- res$ce_result %>%
    pull(outcomes, intervention)
  
  return(list(costs = costs, outcomes = outcomes))
  
}

# psa_res <- lapply(split(psa_input, 1:nrow(psa_input)), FUN = RunPSA, base_case = model) %>%
#   pmap(., bind_rows)

RunPSA <- function(base_case, psalist) {
  
  # run PSA: split PSA inputs into lists by row, calculate PSA (map function), and pull costs and effects together into separate tibbles
  psaTraces <- psalist %>%
    transpose() %>% # convert each row of data into a separate list
    map(.f = genPSA, model = base_case) %>% # apply function to each dataframe in the list to get CE pairs
    pmap(., bind_rows) # pull together costs and outcomes into separate tables for WTP function
  
  return(list(costs = psaTraces$costs, outcomes = psaTraces$outcomes))
  
}

# function takes the maxwtp, a vector or dataframe of costs and outcomes for all interventions
# and returns probability that a given intervention is cost effective 
genProbCE <- function(wtp, costs, outcomes) {
  NMB <- outcomes*wtp - costs                   # NMB = DALYs * wtp - COSTS
  maxNMB <- apply(NMB, MARGIN = 1, max)         # find maximum for each row (nax nmb)
  whichMaxNMB <- (NMB == maxNMB)                # check which intervention has max NMB at that WTP, it is the most cost effective at that WTP
  return(apply(whichMaxNMB, 2, mean))           # calculate probability of cost-effectiveness
}


RunCEAC <- function(costs, outcomes, wtp, by) {
  
  # probability of cost effectiveness calculations
  probCE <- map_dfr(seq(0, wtp, by = by), .f = genProbCE, costs = costs, outcomes = outcomes) %>%
    bind_rows() %>%
    bind_cols(WTP = seq(0, wtp, by = by), .)
  
  return(probCE)
  
}



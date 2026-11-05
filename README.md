# Schizo-Model

Cost-effectiveness of five antipsychotics for schizophrenia in Uganda: chlorpromazine, haloperidol, quetiapine, risperidone and olanzapine.

This is the model behind:

> Lubinga SJ, Mutamba BB, Nganizi A, Babigumira JB. A cost-effectiveness analysis of antipsychotics for treatment of schizophrenia in Uganda. *Applied Health Economics and Health Policy*. 2015;13(5):493–506. [doi:10.1007/s40258-015-0176-3](https://doi.org/10.1007/s40258-015-0176-3)

**Live app:** https://theejonslo.shinyapps.io/schizo-model/

It is a 10-state Markov cohort model with 70 annual cycles. Outcomes are costs (2012 US$) and disability-adjusted life years (DALYs) per patient.

## Model structure

![Structure of the Markov model: residual phase, stopped and relapse states for the first- and second-line drugs, the third-line drug, and death](www/model_structure.svg)

Each year people stay in their state (the loops) or move along one arrow. Everyone starts in the residual phase on their first-line drug. Side effects or lack of efficacy move people to the next drug, as does remission after a relapse on treatment; remission after a relapse off treatment restarts the same drug. The third-line drug is fluphenazine decanoate. The relapse states are spent in hospital.

## How costs are calculated

Each year (cycle length $\tau$ = 1 year) a person's cost depends on the state they are in. Costs are discounted and summed over the cohort trace $\pi_{t,s}$, the share of the cohort in state $s$ in cycle $t$:

$$\text{Cost} = \sum_{t=0}^{T-1} \frac{1}{(1+r)^{t}} \sum_{s} \pi_{t,s}\, C_s$$

**Side effects.** On drug $j$, the probability of extrapyramidal symptoms (EPS) comes from the odds ratio against placebo, and the probability of diabetes from the relative risk against haloperidol:

$$\text{EPS}_j = \frac{\text{OR}_j\, p_0}{1 - p_0 + \text{OR}_j\, p_0} \qquad \text{DM}_j = p_{\text{DM}}\, \text{RR}_j$$

**Staff pay.** Each staff member's hourly pay is their annual salary over the hours worked in a year, $h_k = S_k / H$, for psychiatrists (psy), psychiatric clinical officers (co), psychiatric nurses (nu) and pharmacy assistants (ph).

**A cycle in hospital** (the four relapse states; a stay of $L$ days on chlorpromazine). Ward rounds are counted per week, so they are divided by 7 for a daily cost:

$$C^{\text{IP}}_{\text{staff}} = \left[\frac{1}{7}\left(R_{\text{maj}} \frac{d_{\text{maj}}}{60} \big(N^{\text{maj}}_{\text{psy}} h_{\text{psy}} + N^{\text{maj}}_{\text{co}} h_{\text{co}} + N^{\text{maj}}_{\text{nu}} h_{\text{nu}}\big) + R_{\text{min}} \frac{d_{\text{min}}}{60} \big(N^{\text{min}}_{\text{co}} h_{\text{co}} + N^{\text{min}}_{\text{nu}} h_{\text{nu}}\big)\right) + N_{\text{adm}} \frac{d_{\text{adm}}}{60} h_{\text{nu}}\right] \tau L$$

$$C^{\text{IP}}_{\text{lab}} = \frac{\text{CPI}_{2012}}{\text{CPI}_{2001}} \Big(c_{\text{FBC}} f_{\text{FBC}} + c_{\text{ESR}} f_{\text{ESR}} + p_{\text{TPHA}}\, c_{\text{TPHA}} f_{\text{TPHA}} + (1 - p_{\text{TPHA}})\, c_{\text{VDRL}} f_{\text{VDRL}} + c_{\text{HIV}} f_{\text{HIV}}\Big)\, 12\tau$$

$$C^{\text{IP}}_{\text{capital}} = k_{\text{IP}}\, \tau L \qquad C^{\text{IP}}_{\text{travel}} = 2 c_{\text{trip}}\, p_{\text{carer}} \frac{\tau L\, v_{\text{carer}}}{7} + 2 c_{\text{trip}} \qquad C^{\text{IP}}_{\text{work}} = \lambda (1 - u)\, \text{GDP}\, \tau$$

$$C^{\text{acute}} = L\, c_{\text{CPZ}} n_{\text{CPZ}} + \text{EPS}_{\text{CPZ}}\, L\, c_{\text{bz}} n_{\text{bz}} + \text{DM}_{\text{CPZ}}\, C_{\text{DM}} + C^{\text{IP}}_{\text{staff}} + C^{\text{IP}}_{\text{lab}} + C^{\text{IP}}_{\text{capital}} + C^{\text{IP}}_{\text{travel}} + C^{\text{IP}}_{\text{work}} + C_{\text{carer}}$$

**A cycle as an outpatient** ($V$ visits a month, so $12V\tau$ visits a cycle). The consultation is with a psychiatrist (a share $\sigma_{\text{psy}}$ of patients) or a clinical officer, and a nurse; then the drugs are dispensed:

$$C^{\text{OP}}_{\text{staff}} = \left[\frac{d_{\text{con}}}{60} \big(N_{\text{psy}}\, \sigma_{\text{psy}} h_{\text{psy}} + N_{\text{co}} (1 - \sigma_{\text{psy}}) h_{\text{co}} + N_{\text{nu}}\, \sigma_{\text{nu}} h_{\text{nu}}\big) + N_{\text{ph}} \frac{d_{\text{disp}}}{60} h_{\text{ph}}\right] 12 V \tau$$

$$C^{\text{OP}}_{\text{lab}} = \frac{\text{CPI}_{2012}}{\text{CPI}_{2001}} \big(c_{\text{FBC}} f_{\text{FBC}} + c_{\text{ESR}} f_{\text{ESR}}\big)\, 12\tau \qquad C^{\text{OP}}_{\text{capital}} = k_{\text{OP}}\, 12 V \tau$$

$$C^{\text{OP}}_{\text{travel}} = 2 c_{\text{trip}} (p_{\text{carer}} + 1)\, 12 V \tau \qquad C^{\text{OP}}_{\text{work}} = \lambda (1 - u) \frac{\text{GDP}}{264}\, 12 V \tau$$

On drug $j$ (tablets a day $n_j$ at price $c_j$; benzhexol for EPS for $\delta_{\text{EPS}}$ years), and on fluphenazine decanoate (injections a month):

$$C^{\text{res}}_j = 365\, c_j n_j + \text{EPS}_j\, 365\, c_{\text{bz}} n_{\text{bz}}\, \delta_{\text{EPS}} + \text{DM}_j\, C_{\text{DM}} + C^{\text{OP}}_{\text{staff}} + C^{\text{OP}}_{\text{lab}} + C^{\text{OP}}_{\text{capital}} + C^{\text{OP}}_{\text{travel}} + C^{\text{OP}}_{\text{work}} + C_{\text{carer}}$$

$$C^{\text{res}}_{\text{FLU}} = 12\, c_{\text{FLU}} n_{\text{FLU}} + \text{EPS}_{\text{FLU}}\, 365\, c_{\text{bz}} n_{\text{bz}}\, \delta_{\text{EPS}} + p_{\text{DM,FLU}}\, C_{\text{DM}} + C^{\text{OP}}_{\text{staff}} + C^{\text{OP}}_{\text{lab}} + C^{\text{OP}}_{\text{capital}} + C^{\text{OP}}_{\text{travel}} + C^{\text{OP}}_{\text{work}} + C_{\text{carer}}$$

**Carers** lose a share $\ell$ of their work time (264 work days a year), in hospital and as an outpatient: $C_{\text{carer}} = \lambda (1 - u)\, \text{GDP}\, \tau\, \ell\, \tfrac{264}{365}$. Time off treatment (the stopped states) costs nothing.

| Symbol | Input in `data_input.xlsx` | Base case |
|---|---|---|
| $S_{\text{psy}}, S_{\text{co}}, S_{\text{nu}}, S_{\text{ph}}$; $H$ | `rate_psy`, `rate_pco`, `rate_pno`, `rate_pharm`; `workhours_year` | $5,049, $1,261, $970, $1,087; 2,000 hours |
| $R_{\text{maj}}, R_{\text{min}}$; $d_{\text{maj}}, d_{\text{min}}$ | `major_week`, `minor_week`; `duration_major`, `duration_minor` | 1 and 2 rounds a week; 20 minutes each |
| $N^{\text{maj}}_{\text{psy}}, N^{\text{maj}}_{\text{co}}, N^{\text{maj}}_{\text{nu}}$; $N^{\text{min}}_{\text{co}}, N^{\text{min}}_{\text{nu}}$ | `PSY_major`, `PCO_major`, `PNO_major`; `PCO_minor`, `PNO_minor` | 1, 1, 2; 1, 2 |
| $N_{\text{adm}}$, $d_{\text{adm}}$ | `PNO_drugadmin`, `duration_drugadmin` | 1 nurse, 20 minutes a day |
| $N_{\text{psy}}, N_{\text{co}}, N_{\text{nu}}, N_{\text{ph}}$; $\sigma_{\text{psy}}, \sigma_{\text{nu}}$ | `PSY`, `PCO`, `nurse`, `pharm`; `PSY_prop`, `nurse_prop` | 1 each; 0.25, 1 |
| $d_{\text{con}}$, $d_{\text{disp}}$; $V$ | `duration_consult`, `duration_drugdisp`; `visitation_frequency` | 20 and 10 minutes; 1 visit a month |
| $L$; $\tau$; $r$ | `LOS`; `time`; `discCosts` | 60 days; 1 year; 3% |
| $c$ and $f$ for each test; $p_{\text{TPHA}}$ | `FBC_cost`, `FBC_freq`, … `HIV_cost`, `HIV_freq`; `TPHA_proportion` | tests once a month; 0.1 |
| $k_{\text{IP}}$, $k_{\text{OP}}$ | `capital_costs_per_inpatient_per_day`, `capital_costs_per_outpatient_visit` | $8.46 a day, $2.49 a visit |
| $c_{\text{trip}}$; $p_{\text{carer}}$; $v_{\text{carer}}$ | `averagetravel_cost`; `caregiver_prop`; `caregiver_freq` (per `every` = 7 days) | $2 each way; 0.4; 1 visit a week |
| $\lambda$, $u$, GDP; $\ell$ | `laborforcerate`, `unemploymentrate`, `GDP_percapita`; `DALY_caretaker` | 0.706, 0.052, $547; 0.25 |
| $p_0$; $p_{\text{DM}}$, $p_{\text{DM,FLU}}$; $C_{\text{DM}}$ | `pEPSPLA`; `pDiabetesHAL`, `pDiabetesFLU`; `cDiabetes` | 0.064; 0.052, 0.06; $627 a year |
| $\text{OR}_j$, $\text{RR}_j$; $c_j$, $n_j$ | `eps…`, `diab…`; `…_cost`, `…_qty` for each drug | Table 2 of the paper |
| $c_{\text{bz}}$, $n_{\text{bz}}$; $\delta_{\text{EPS}}$ | `artane_cost`, `artane_qty`; `durEPS` | $0.016, 1 a day; 0.25 years |

The functions `ip_costs()`, `op_costs()`, `residual_cost()`, `acute_cost()` and `third_line_cost()` in `R/functions.R` follow these equations term by term.

## Contents

| Path | What it is |
|---|---|
| `new_model_10_21_15.xlsm` | The original Excel model (2015) |
| `data/data_input.xlsx` | All inputs: parameters, transition probabilities, CPI, life tables |
| `R/functions.R` | Model, cost-effectiveness frontier, one-way and probabilistic sensitivity analyses |
| `R/data_input.R` | Reads the inputs and model settings into one list |
| `R/00_packages.R` | Loads the packages the model uses |
| `main.R` | Runs the base case, one-way analysis (tornado), PSA and acceptability curves |
| `app.R`, `app/`, `www/` | The interactive Shiny app: page layout and server (`app.R`), helpers, charts and input editors (`app/`), styles, font, icon and the model structure figure (`www/`) |

## Changes since the 2015 publication

The published results came from the Excel model, which is kept in this repository unchanged. The R version reproduces its structure, with the corrections below, so its results differ from the paper.

| Change | Excel model (2015) | R model | Effect on base-case costs |
|---|---|---|---|
| **Cost of acute (hospitalised) states** | The four acute-state cost cells multiplied by the YLD columns (state occupancy × disability weight) instead of the state occupancy, so hospital costs were scaled down by the disability weight (about 0.77) | Costs are multiplied by state occupancy | About +$310 to +$415 per patient |
| **Laboratory-cost inflation** | Lab costs were inflated with a mistyped 2012 CPI (196.64 instead of 196.46) | All CPI adjustments are read from the CPI sheet of `data_input.xlsx` | About −$0.65 per patient |
| **HIV test in inpatient lab costs** | The HIV-test frequency multiplied the VDRL term, so the HIV test was never costed | HIV test cost × frequency is added as its own term | About +$50 to +$68 per patient |
| **Drug price years** | Chlorpromazine and haloperidol prices (MSH International Drug Price Indicator Guide 2013, already in US$) were inflated with an index of 145.14, which is not in the CPI table; inpatient chlorpromazine was inflated from 2010 instead. Benzhexol (Artane), already converted to 2012 US$ from its NMS October 2010 price, was inflated a second time with 196.46/145.14 | MSH prices are used as published, like olanzapine and risperidone. Benzhexol is converted once: 15,750 UGX per 500 tablets, at 2,577 UGX/$, inflated from 2010 | −$6 to −$17 |
| **Haloperidol dose** | One 10 mg tablet a day at $0.0206 | Two 5 mg tablets a day at $0.022, as in Table 2 of the paper | +$9 to +$12 (haloperidol, and the three drugs that use it second line) |
| **Quetiapine price** | $0.388 per 200 mg tablet | $0.188, the pharmacy-survey price in Table 2 of the paper. The source of $0.388 could not be traced | −$276 (quetiapine only) |
| **Benzhexol in residual states** (bug) | Outpatients with extrapyramidal symptoms (EPS) were costed for the 60-day inpatient supply of benzhexol. The cell for an outpatient supply was never used | Costed for the duration of EPS (`durEPS`, 0.25 years), the duration used for the disability from EPS | +$0.4 to +$0.5 |
| **Background mortality** (bug) | The WHO life table for Uganda gives the probability of dying within each 5-year age group (₅qₓ). These were used as yearly probabilities, so people died about five times too fast: life expectancy from age 25 was 16.7 years, against 39 years in the same life table | The 5-year probabilities are converted to yearly ones, 1 − (1 − ₅qₓ)^(1/5), using the width of each age group. Life expectancy from age 25 is now 33.3 years, including the extra mortality from schizophrenia | +$2,619 to +$2,657; DALYs −5.54 to −5.55 |
| **Staff time** (bug) | Ward rounds are counted per week, but only the psychiatrist's major round was divided by 7 to give a daily cost, and the daily drug round was added in hours rather than dollars. In outpatient care, the consultation length applied to the psychiatrist only, so the clinical officer and nurse were costed for a full hour | Every weekly round is divided by 7, the drug round is costed at the nurse's hourly pay, and the consultation length applies to everyone at it (equations above). Staff time per hospital stay falls from $123.25 to $30.63, and per year of outpatient care from $15.11 to $7.44 | −$257 to −$292 |
| **Probability of cost-effectiveness** (bug) | Net benefit was computed as WTP × DALYs − cost, with DALYs incurred, so more DALYs counted as a benefit (the PSA sheet, and the net-benefit equation in the paper). The acceptability curve and the published shares (risperidone most cost-effective in 85% of simulations at $547 and 74% at $1,641) used this | DALYs count as a loss: net benefit = −WTP × DALYs − cost | None on the base case; changes the acceptability curve |

Inpatient drug costs follow the length of stay, as in the Excel model; this matters only in the sensitivity analyses. The R model runs 70 yearly cycles (to age 94), the Excel model 62 (to age 86); with the corrected mortality about 2% of the cohort is still alive at 86, so this adds about $4 and 0.02 DALYs per patient. The fluphenazine decanoate price ($1.406 per 25 mg ampoule, from the NMS October 2010 catalogue) is kept as published, although it cannot be reproduced from the catalogue price of 3,950 UGX per ampoule at 2,577 UGX/$. The schizophrenia-related death probabilities (WHO 2008 deaths by cause, in wider age bands) are used as yearly probabilities, as before.

Base-case results per patient (costs in 2012 US$):

| Drug | Cost, Excel model (2015) | Cost, R model | DALYs, Excel model | DALYs, R model |
|---|---|---|---|---|
| Chlorpromazine | $4,617 | $7,449 | 27.608 | 22.066 |
| Haloperidol | $4,440 | $7,195 | 27.563 | 22.012 |
| Quetiapine | $4,943 | $7,441 | 27.562 | 22.009 |
| Risperidone | $4,424 | $7,178 | 27.557 | 22.008 |
| Olanzapine | $4,454 | $7,188 | 27.552 | 22.003 |

In the Excel model, olanzapine cost $5,868 per DALY averted compared with risperidone. In the R model, chlorpromazine, haloperidol and quetiapine are dominated. Risperidone costs the least; olanzapine costs $9.79 more per patient and averts 0.005 DALYs, or $1,999 per DALY averted. That is above three times Uganda's 2012 GDP per capita ($1,641), so risperidone is the cost-effective choice at that threshold. The differences are small, though: in the probabilistic analysis (1,000 draws), risperidone is the most cost-effective drug in about 55% of draws at $1,641 per DALY averted, and olanzapine in about 44%.

## Running it

Install the packages listed at the top of `R/00_packages.R`, then open `Schizo-Model.Rproj` in RStudio and run `main.R`.

The model functions take everything they need as arguments:

```r
inputs <- read_inputs()   # input tables, plus nStates, nCycles and interventions
base_case <- RunModel(inputs,
                      baseparms = getparms(inputs$parameters, value = basecase),
                      basetransitions = getparms(inputs$transitions, value = basecase))
psa_res <- RunPSA(base_case, nsims = 1000)
```

The result of `RunModel()` carries its inputs, so `owsa()` and `RunPSA()` only need the model object.

## The app

The Shiny app runs the same model. You can:

- edit any input: parameters and their distributions, transition probabilities, mortality and life-expectancy tables, the price index, the time horizon and the willingness-to-pay threshold, then download your inputs as a spreadsheet and upload them again later;
- see the base case, with the efficiency frontier on a cost-effectiveness plane;
- run the one-way analysis for any pair of drugs, as a tornado diagram of the incremental cost, DALYs averted, cost per DALY averted or net monetary benefit;
- run the probabilistic analysis, with acceptability curves for all five drugs and a cost-effectiveness plane for each drug against a comparator;
- read how this version differs from the published model.

To run it locally, install the packages the model uses (listed in `R/00_packages.R`) plus the app's own:

```r
install.packages(c("shiny", "bslib", "reactable", "writexl"))
shiny::runApp()
```

Shiny loads the files in `R/` automatically before `app.R`. To deploy it to shinyapps.io from this folder:

```r
rsconnect::deployApp(appName = "schizo-model")
```

`.rscignore` leaves the Excel model out of the deployed bundle. The app uses the Atkinson Hyperlegible Next and Mono fonts (SIL Open Font License, `www/fonts/OFL.txt`), bundled so they work offline.

## Licence

The code is released under the MIT License (see `LICENSE`). The bundled fonts are under the SIL Open Font License (see `www/fonts/OFL.txt`). The published article is not covered by this licence.

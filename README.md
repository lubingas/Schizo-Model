# Schizo-Model

Cost-effectiveness of five antipsychotics for schizophrenia in Uganda: chlorpromazine, haloperidol, quetiapine, risperidone and olanzapine.

This is the model behind:

> Lubinga SJ, Mutamba BB, Nganizi A, Babigumira JB. A cost-effectiveness analysis of antipsychotics for treatment of schizophrenia in Uganda. *Applied Health Economics and Health Policy*. 2015;13(5):493–506. [doi:10.1007/s40258-015-0176-3](https://doi.org/10.1007/s40258-015-0176-3)

**Live app:** https://theejonslo.shinyapps.io/schizo-model/

It is a 10-state Markov cohort model with 70 annual cycles. Outcomes are costs (2012 US$) and disability-adjusted life years (DALYs) per patient.

## Model structure

![Structure of the Markov model: residual phase, stopped and relapse states for the first- and second-line drugs, the third-line drug, and death](www/model_structure.svg)

Each year people stay in their state (the loops) or move along one arrow. Everyone starts in the residual phase on their first-line drug. Side effects or lack of efficacy move people to the next drug, as does remission after a relapse on treatment; remission after a relapse off treatment restarts the same drug. The third-line drug is fluphenazine decanoate. The relapse states are spent in hospital.

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
| **Probability of cost-effectiveness** (bug) | Net benefit was computed as WTP × DALYs − cost, with DALYs incurred, so more DALYs counted as a benefit (the PSA sheet, and the net-benefit equation in the paper). The acceptability curve and the published shares (risperidone most cost-effective in 85% of simulations at $547 and 74% at $1,641) used this | DALYs count as a loss: net benefit = −WTP × DALYs − cost | None on the base case; changes the acceptability curve |

Inpatient drug costs follow the length of stay, as in the Excel model; this matters only in the sensitivity analyses. The R model runs 70 yearly cycles (to age 94), the Excel model 62 (to age 86); with the corrected mortality about 2% of the cohort is still alive at 86, so this adds about $4 and 0.02 DALYs per patient. The fluphenazine decanoate price ($1.406 per 25 mg ampoule, from the NMS October 2010 catalogue) is kept as published, although it cannot be reproduced from the catalogue price of 3,950 UGX per ampoule at 2,577 UGX/$. The schizophrenia-related death probabilities (WHO 2008 deaths by cause, in wider age bands) are used as yearly probabilities, as before.

Base-case results per patient (costs in 2012 US$):

| Drug | Cost, Excel model (2015) | Cost, R model | DALYs, Excel model | DALYs, R model |
|---|---|---|---|---|
| Chlorpromazine | $4,617 | $7,741 | 27.608 | 22.066 |
| Haloperidol | $4,440 | $7,462 | 27.563 | 22.012 |
| Quetiapine | $4,943 | $7,708 | 27.562 | 22.009 |
| Risperidone | $4,424 | $7,443 | 27.557 | 22.008 |
| Olanzapine | $4,454 | $7,445 | 27.552 | 22.003 |

In the Excel model, olanzapine cost $5,868 per DALY averted compared with risperidone. In the R model, chlorpromazine, haloperidol and quetiapine are dominated. Risperidone costs the least; olanzapine costs $1.63 more per patient and averts 0.005 DALYs, or $332 per DALY averted. That is below three times Uganda's 2012 GDP per capita ($1,641), so olanzapine is the cost-effective choice at that threshold, but the differences are small enough that the probabilistic analysis matters more than the point estimate.

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

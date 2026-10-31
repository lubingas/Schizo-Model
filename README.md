# Schizo-Model

Cost-effectiveness of five antipsychotics for schizophrenia in Uganda: chlorpromazine, haloperidol, quetiapine, risperidone and olanzapine.

This is the model behind:

> Lubinga SJ, Mutamba BB, Nganizi A, Babigumira JB. A cost-effectiveness analysis of antipsychotics for treatment of schizophrenia in Uganda. *Applied Health Economics and Health Policy*. 2015;13(5):493–506. [doi:10.1007/s40258-015-0176-3](https://doi.org/10.1007/s40258-015-0176-3)

It is a 10-state Markov cohort model with 70 annual cycles. Outcomes are costs (2012 US$) and disability-adjusted life years (DALYs) per patient.

## Contents

| Path | What it is |
|---|---|
| `new_model_10_21_15.xlsm` | The original Excel model (2015) |
| `data/data_input.xlsx` | All inputs: parameters, transition probabilities, CPI, life tables |
| `R/functions.R` | Model, cost-effectiveness frontier, one-way and probabilistic sensitivity analyses |
| `R/data_input.R` | Reads the inputs |
| `main.R` | Runs the base case, one-way analysis (tornado), PSA and acceptability curves |

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

Inpatient drug costs follow the length of stay, as in the Excel model; this matters only in the sensitivity analyses. The fluphenazine decanoate price ($1.406 per 25 mg ampoule, from the NMS October 2010 catalogue) is kept as published, although it cannot be reproduced from the catalogue price of 3,950 UGX per ampoule at 2,577 UGX/$.

Base-case results, cost per patient (2012 US$):

| Drug | Excel model (2015) | R model | DALYs |
|---|---|---|---|
| Chlorpromazine | $4,617 | $5,084 | 27.608 |
| Haloperidol | $4,440 | $4,843 | 27.563 |
| Quetiapine | $4,943 | $5,070 | 27.562 |
| Risperidone | $4,424 | $4,817 | 27.557 |
| Olanzapine | $4,454 | $4,816 | 27.552 |

In the Excel model, olanzapine cost $5,868 per DALY averted compared with risperidone. In the R model, olanzapine has the lowest cost and the fewest DALYs, so it dominates the other four drugs. Its cost advantage over risperidone is under $1 per patient, so the probabilistic analysis matters more than the point estimate.

## Running it

Open `Schizo-Model.Rproj` in RStudio and run `main.R`.

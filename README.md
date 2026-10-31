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

The R version reproduces the Excel model's structure, with these corrections. Results therefore differ from the published paper.

| Change | Before | After | Effect on base-case costs |
|---|---|---|---|
| **Laboratory-cost inflation** | Lab costs were inflated with a mistyped 2012 CPI (196.64 instead of 196.46) | CPI adjustments are read from the CPI sheet of `data_input.xlsx` (except Artane, below) | About −$0.65 per patient |
| **HIV test in inpatient lab costs** | The HIV-test frequency multiplied the VDRL term, so the HIV test was never costed | HIV test cost × frequency is added as its own term | About +$50 to +$68 per patient |
| **Cost-effectiveness frontier** | Options were walked in order of effectiveness, and dominated options could stay on the frontier | Standard method: remove strongly dominated options, then extendedly dominated ones, then compute ICERs along the frontier | Changes which options are reported as cost-effective |

With these changes, olanzapine has the lowest cost and the fewest DALYs in the base case, so it dominates the other four drugs. Its cost advantage over risperidone is about $1 per patient, so the probabilistic analysis matters more than the point estimate.

One open item: the Artane (benzhexol) price is inflated with an index of 145.14, which is not in the CPI table (2009 = 144.24, 2010 = 149.97). It is kept as in the Excel model until the source year is confirmed.

## Running it

Open `Schizo-Model.Rproj` in RStudio and run `main.R`.

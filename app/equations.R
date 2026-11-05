# app/equations.R
# ---------------------------
# How the model calculates results, for the About tab, typeset with MathJax.
# The same equations, with a table of the inputs behind each symbol, are in the README;
# mtrace(), RunModel(), state_weight(), third_line_weight(), ip_costs(), op_costs(), residual_cost(),
# acute_cost(), third_line_cost(), IncrementalAnalysis(), genProbCE(), owsaPair(), getTPbounds(), owsa()
# and draw_psa_df() in R/functions.R follow them.
# ---------------------------

eq <- function(tex) div(class = "eq", paste0("\\[", tex, "\\]"))

model_equations <- function() {
  withMathJax(
    p("The model follows a cohort of patients aged \\(a_0\\) = 25 for \\(T\\) = 70 yearly cycles (cycle length \\(\\tau\\) = 1 year), so in cycle ",
      "\\(t = 0, 1, \\dots, T-1\\) they are aged \\(a_t = a_0 + t\\tau\\). Each of the five drugs is the first-line drug \\(j_1\\) of a treatment sequence. ",
      "The second-line drug \\(j_2\\) is risperidone after chlorpromazine or haloperidol, and haloperidol after quetiapine, risperidone or olanzapine. ",
      "Fluphenazine decanoate (FLU) is the third line for every patient."),

    h3("States and transitions"),
    tags$table(class = "diff-table",
      tags$thead(tags$tr(tags$th("State"), tags$th("Description"))),
      tags$tbody(
        tags$tr(tags$td("\\(\\text{R}_1, \\text{R}_2, \\text{R}_3\\)"), tags$td("Residual phase on the first-line drug, the second-line drug, or fluphenazine decanoate")),
        tags$tr(tags$td("\\(\\text{S}_1, \\text{S}_2\\)"), tags$td("Residual phase after stopping the first- or second-line drug (untreated)")),
        tags$tr(tags$td("\\(\\text{A}'_1, \\text{A}'_2\\)"), tags$td("Relapse after stopping the first- or second-line drug")),
        tags$tr(tags$td("\\(\\text{A}_1, \\text{A}_2\\)"), tags$td("Relapse while on the first- or second-line drug")),
        tags$tr(tags$td("\\(\\text{D}\\)"), tags$td("Dead")))),
    p("Each year, a patient in the residual phase on drug \\(j\\) relapses with probability \\(x^{\\text{rel}}_j\\), switches to the next drug because of ",
      "side effects (\\(x^{\\text{AE}}_j\\)) or lack of efficacy (\\(x^{\\text{eff}}_j\\)), stops for other reasons (\\(x^{\\text{oth}}_j\\)), or continues ",
      "(\\(x^{\\text{cont}}_j\\)). For each drug these five come from the same trial counts and add up to 1. Patients who have stopped treatment relapse with ",
      "probability \\(x^{\\text{rel}}_0\\). Relapses are treated in hospital and remit with probability \\(\\rho\\): after a relapse off treatment the patient ",
      "restarts the same drug, and after a relapse on treatment moves to the next one. Every living patient dies with probability \\(q_t\\) in cycle \\(t\\), ",
      "whatever their state, and the rest stay where they are. For lines \\(\\ell\\) = 1 and 2:"),
    eq(r"(P_t(\text{R}_\ell, \text{R}_{\ell+1}) = x^{\text{AE}}_{j_\ell} + x^{\text{eff}}_{j_\ell} \qquad P_t(\text{R}_\ell, \text{S}_\ell) = x^{\text{oth}}_{j_\ell} \qquad P_t(\text{R}_\ell, \text{A}_\ell) = x^{\text{rel}}_{j_\ell})"),
    eq(r"(P_t(\text{S}_\ell, \text{A}'_\ell) = x^{\text{rel}}_0 \qquad P_t(\text{A}'_\ell, \text{R}_\ell) = \rho \qquad P_t(\text{A}_\ell, \text{R}_{\ell+1}) = \rho)"),
    eq(r"(P_t(s, \text{D}) = q_t \ \ (s \ne \text{D}) \qquad P_t(s, s) = 1 - \sum_{s' \ne s} P_t(s, s') \qquad P_t(\text{D}, \text{D}) = 1)"),
    p("So a patient stays on drug \\(j_\\ell\\) with probability \\(x^{\\text{cont}}_{j_\\ell} - q_t\\), and stays on fluphenazine decanoate until death. ",
      "Every patient starts in \\(\\text{R}_1\\), and the cohort trace \\(\\pi_t\\), the share of the cohort in each state, moves one cycle at a time:"),
    eq(r"(\pi_0 = (1, 0, \dots, 0) \qquad \pi_{t+1} = \pi_t\, P_t)"),

    h3("Deaths"),
    p("The WHO life table for Uganda gives \\({}_5q_x\\), the probability of dying within the 5-year age group that starts at age \\(x\\), which is converted ",
      "to a yearly probability. Deaths from schizophrenia, \\(q^{\\text{SZ}}_x\\) (WHO 2008 deaths by cause, in wider age groups), are used as yearly probabilities. ",
      "A patient dies of one cause or the other, with each table read at the age group that contains \\(a_t\\):"),
    eq(r"(q^{\text{pop}}_x = 1 - (1 - {}_5q_x)^{1/5} \qquad q_t = q^{\text{pop}}_{a_t} + q^{\text{SZ}}_{a_t} - q^{\text{pop}}_{a_t}\, q^{\text{SZ}}_{a_t})"),

    h3("Side effects"),
    p("On drug \\(j\\), extrapyramidal symptoms (EPS) come from the odds ratio against placebo, and diabetes from the relative risk against haloperidol:"),
    eq(r"(\text{EPS}_j = \frac{\text{OR}_j\, p_0}{1 - p_0 + \text{OR}_j\, p_0} \qquad \text{DM}_j = p_{\text{DM}}\, \text{RR}_j)"),

    h3("DALYs"),
    p("Each year, the patients who died in the previous cycle lose the life expectancy \\(e(a_t)\\) at their age (from the same life table), and the living ",
      "accrue the disability weight \\(w_s\\) of their state. Both are discounted at rate \\(r_{\\text{O}}\\):"),
    eq(r"(\text{YLL}_t = q_{t-1}\, e(a_t) \sum_{s \ne \text{D}} \pi_{t-1,s} \ \ (\text{YLL}_0 = 0) \qquad \text{YLD}_t = \sum_s \pi_{t,s}\, w_s \qquad \text{DALYs} = \sum_{t=0}^{T-1} \frac{\text{YLL}_t + \text{YLD}_t}{(1 + r_{\text{O}})^{t}})"),
    p("A state's disability weight is that of its phase (residual or acute) plus the expected disability from side effects: weight gain (with probability ",
      "\\(p^{\\text{wt}}_j\\)), EPS (for \\(\\delta_{\\text{EPS}}\\) years) and diabetes. Patients who have stopped treatment have no side effects, and in hospital ",
      "every patient is on chlorpromazine (\\(w_{\\text{A}}\\) is the weight of all four relapse states):"),
    eq(r"(w_{\text{R}_\ell} = W_{\text{res}} + p^{\text{wt}}_{j_\ell} W_{\text{wt}} + \text{EPS}_{j_\ell}\, \delta_{\text{EPS}}\, W_{\text{EPS}} + \text{DM}_{j_\ell} W_{\text{DM}} \qquad w_{\text{R}_3} = W_{\text{res}} + p^{\text{wt}}_{\text{FLU}} W_{\text{wt}} + \text{EPS}_{\text{FLU}}\, \delta_{\text{EPS}}\, W_{\text{EPS}} + p_{\text{DM,FLU}} W_{\text{DM}})"),
    eq(r"(w_{\text{S}_\ell} = W_{\text{res}} \qquad w_{\text{A}} = W_{\text{acute}} + p^{\text{wt}}_{\text{CPZ}} W_{\text{wt}} + \text{EPS}_{\text{CPZ}}\, \delta_{\text{EPS}}\, W_{\text{EPS}} + \text{DM}_{\text{CPZ}} W_{\text{DM}} \qquad w_{\text{D}} = 0)"),
    p("The disability from diabetes includes its complications:"),
    eq(r"(W_{\text{DM}} = W_{\text{diab}} + p_{\text{neuro}} W_{\text{neuro}} + p_{\text{foot}} W_{\text{foot}} + p_{\text{imp}} W_{\text{imp}} + p_{\text{neph}} W_{\text{neph}})"),
    p("Every living patient has the same probability of dying, so years of life lost are the same for all five drugs, and the drugs differ in DALYs only through disability."),

    h3("Costs"),
    p("Each year a patient's cost depends on their state: \\(C^{\\text{res}}_{j_\\ell}\\) in \\(\\text{R}_\\ell\\), \\(C^{\\text{res}}_{\\text{FLU}}\\) in \\(\\text{R}_3\\), ",
      "\\(C^{\\text{acute}}\\) in the four relapse states, and nothing in the stopped states or once dead. Costs are discounted at rate \\(r_{\\text{C}}\\):"),
    eq(r"(\text{Cost} = \sum_{t=0}^{T-1} \frac{1}{(1+r_{\text{C}})^{t}} \sum_{s} \pi_{t,s}\, C_s)"),
    p(tags$strong("Staff pay."), " Hourly pay is annual salary over the hours worked in a year, \\(h_k = S_k / H\\), for psychiatrists (psy), psychiatric ",
      "clinical officers (co), psychiatric nurses (nu) and pharmacy assistants (ph)."),
    p(tags$strong("A cycle in hospital."), " The four relapse states: a stay of \\(L\\) days on chlorpromazine. Ward rounds are counted per week, so they are ",
      "divided by 7 for a daily cost. Laboratory tests (a blood count, syphilis and HIV) are done once per stay. A share ",
      "\\(p^{\\text{IP}}_{\\text{carer}}\\) of patients have a carer visit \\(v_{\\text{carer}}\\) times a week."),
    eq(r"(C^{\text{IP}}_{\text{staff}} = \left[\frac{1}{7}\left(R_{\text{maj}} \frac{d_{\text{maj}}}{60} \big(N^{\text{maj}}_{\text{psy}} h_{\text{psy}} + N^{\text{maj}}_{\text{co}} h_{\text{co}} + N^{\text{maj}}_{\text{nu}} h_{\text{nu}}\big) + R_{\text{min}} \frac{d_{\text{min}}}{60} \big(N^{\text{min}}_{\text{co}} h_{\text{co}} + N^{\text{min}}_{\text{nu}} h_{\text{nu}}\big)\right) + N_{\text{adm}} \frac{d_{\text{adm}}}{60} h_{\text{nu}}\right] \tau L)"),
    eq(r"(C^{\text{IP}}_{\text{lab}} = \frac{\text{CPI}_{2012}}{\text{CPI}_{2001}} \Big(c_{\text{FBC}} f_{\text{FBC}} + c_{\text{ESR}} f_{\text{ESR}} + p_{\text{TPHA}}\, c_{\text{TPHA}} f_{\text{TPHA}} + (1 - p_{\text{TPHA}})\, c_{\text{VDRL}} f_{\text{VDRL}} + c_{\text{HIV}} f_{\text{HIV}}\Big))"),
    eq(r"(C^{\text{IP}}_{\text{capital}} = k_{\text{IP}}\, \tau L \qquad C^{\text{IP}}_{\text{travel}} = 2 c_{\text{trip}}\, p^{\text{IP}}_{\text{carer}} \frac{\tau L\, v_{\text{carer}}}{7} + 2 c_{\text{trip}} \qquad C^{\text{IP}}_{\text{work}} = \lambda (1 - u)\, \text{GDP}\, \tau)"),
    eq(r"(C^{\text{acute}} = L\, c_{\text{CPZ}} n_{\text{CPZ}} + \text{EPS}_{\text{CPZ}}\, L\, c_{\text{bz}} n_{\text{bz}} + \text{DM}_{\text{CPZ}}\, C_{\text{DM}} + C^{\text{IP}}_{\text{staff}} + C^{\text{IP}}_{\text{lab}} + C^{\text{IP}}_{\text{capital}} + C^{\text{IP}}_{\text{travel}} + C^{\text{IP}}_{\text{work}} + C_{\text{carer}})"),
    p(tags$strong("A cycle as an outpatient."), " \\(V\\) visits a month, so \\(12V\\tau\\) visits a cycle. The consultation is with a psychiatrist (a share ",
      "\\(\\sigma_{\\text{psy}}\\) of patients) or a clinical officer, and a nurse; then the drugs are dispensed. Blood counts are done once a cycle. ",
      "A share \\(p^{\\text{OP}}_{\\text{carer}}\\) of patients come to their visits with a carer."),
    eq(r"(C^{\text{OP}}_{\text{staff}} = \left[\frac{d_{\text{con}}}{60} \big(N_{\text{psy}}\, \sigma_{\text{psy}} h_{\text{psy}} + N_{\text{co}} (1 - \sigma_{\text{psy}}) h_{\text{co}} + N_{\text{nu}}\, \sigma_{\text{nu}} h_{\text{nu}}\big) + N_{\text{ph}} \frac{d_{\text{disp}}}{60} h_{\text{ph}}\right] 12 V \tau)"),
    eq(r"(C^{\text{OP}}_{\text{lab}} = \frac{\text{CPI}_{2012}}{\text{CPI}_{2001}} \big(c_{\text{FBC}} f_{\text{FBC}} + c_{\text{ESR}} f_{\text{ESR}}\big) \qquad C^{\text{OP}}_{\text{capital}} = k_{\text{OP}}\, 12 V \tau)"),
    eq(r"(C^{\text{OP}}_{\text{travel}} = 2 c_{\text{trip}} (p^{\text{OP}}_{\text{carer}} + 1)\, 12 V \tau \qquad C^{\text{OP}}_{\text{work}} = \lambda (1 - u) \frac{\text{GDP}}{264}\, 12 V \tau)"),
    p("On drug \\(j\\) (\\(n_j\\) tablets a day at \\(c_j\\) each; benzhexol for EPS for \\(\\delta_{\\text{EPS}}\\) years), and on fluphenazine decanoate (injections a month):"),
    eq(r"(C^{\text{res}}_j = 365\, c_j n_j + \text{EPS}_j\, 365\, c_{\text{bz}} n_{\text{bz}}\, \delta_{\text{EPS}} + \text{DM}_j\, C_{\text{DM}} + C^{\text{OP}}_{\text{staff}} + C^{\text{OP}}_{\text{lab}} + C^{\text{OP}}_{\text{capital}} + C^{\text{OP}}_{\text{travel}} + C^{\text{OP}}_{\text{work}} + C_{\text{carer}})"),
    eq(r"(C^{\text{res}}_{\text{FLU}} = 12\, c_{\text{FLU}} n_{\text{FLU}} + \text{EPS}_{\text{FLU}}\, 365\, c_{\text{bz}} n_{\text{bz}}\, \delta_{\text{EPS}} + p_{\text{DM,FLU}}\, C_{\text{DM}} + C^{\text{OP}}_{\text{staff}} + C^{\text{OP}}_{\text{lab}} + C^{\text{OP}}_{\text{capital}} + C^{\text{OP}}_{\text{travel}} + C^{\text{OP}}_{\text{work}} + C_{\text{carer}})"),
    p(tags$strong("Carers"), " lose a share \\(\\ell\\) of their work time (264 work days a year), in hospital and as an outpatient. ",
      "Time off treatment (the stopped states) costs nothing."),
    eq(r"(C_{\text{carer}} = \lambda (1 - u)\, \text{GDP}\, \tau\, \ell\, \tfrac{264}{365})"),

    h3("Comparing the drugs"),
    p("The drugs are compared on the efficiency frontier. A drug is dominated if another costs no more and incurs no more DALYs, and is better on one of the two. ",
      "The rest are sorted by cost, and a drug is extendedly dominated, and dropped, if its ICER against the drug before it is higher than the ICER of the drug ",
      "after it; this repeats until the ICERs rise along the frontier. Between neighbours on the frontier:"),
    eq(r"(\text{ICER}_k = \frac{\text{Cost}_k - \text{Cost}_{k-1}}{\text{DALYs}_{k-1} - \text{DALYs}_k})"),

    h3("Net benefit and the probability of cost-effectiveness"),
    p("At a willingness to pay of WTP per DALY averted, a drug's net monetary benefit counts its DALYs as a loss. In the probabilistic analysis, a drug's ",
      "probability of being cost-effective is the share of the \\(N\\) draws in which it has the highest net benefit (the acceptability curve plots this against WTP):"),
    eq(r"(\text{NMB}_j = -\text{WTP} \times \text{DALYs}_j - \text{Cost}_j \qquad P_j(\text{WTP}) = \frac{1}{N} \sum_{n=1}^{N} \mathbf{1}\Big[\text{NMB}^{(n)}_j = \max_i \text{NMB}^{(n)}_i\Big])"),
    p("The one-way analysis compares a drug \\(I\\) with a comparator \\(K\\): the incremental cost \\(\\Delta C = \\text{Cost}_I - \\text{Cost}_K\\), the DALYs averted ",
      "\\(\\Delta E = \\text{DALYs}_K - \\text{DALYs}_I\\), the ICER \\(\\Delta C / \\Delta E\\), or the incremental net benefit \\(\\text{WTP} \\times \\Delta E - \\Delta C\\). ",
      "The ICER is shown only when \\(I\\) costs more and averts more; otherwise the run is labelled dominant (\\(\\Delta C \\le 0\\) and \\(\\Delta E \\ge 0\\)), ",
      "dominated (\\(\\Delta C \\ge 0\\) and \\(\\Delta E \\le 0\\)), or cheaper and less effective."),

    h3("Sensitivity analyses"),
    p("In the one-way analysis each parameter is set to its low and its high value. Each transition probability \\(x_i\\) is set to the 2.5% and 97.5% points of ",
      "its own distribution, \\(\\text{Beta}(\\alpha_i, \\sum_k \\alpha_k - \\alpha_i)\\), where \\(\\alpha\\) are the trial counts of that drug's set, and the others in ",
      "the set are rescaled so it still adds up to 1: \\(x_k' = x_k\\, (1 - x_i') / (1 - x_i)\\)."),
    p("In the probabilistic analysis, probabilities, proportions and disability weights are drawn from beta distributions; the odds ratios for EPS and relative ",
      "risks for diabetes from log-normal distributions, with a log-scale standard deviation of \\((\\ln \\text{high} - \\ln \\text{low}) / 3.92\\); costs and prices, ",
      "staff times, salaries and working hours, GDP per capita, the length of stay, visit frequencies and the carers' lost work time from normal distributions; ",
      "and each drug's five transition probabilities from a Dirichlet distribution with its trial counts \\(\\alpha\\). The age, cycle length, discount rates, EPS ",
      "duration, ward-round schedule, staff numbers, doses and test frequencies are fixed."),
    p(class = "small-print", "Every symbol is an input on the Inputs tab; the README on GitHub lists which input each one is.")
  )
}

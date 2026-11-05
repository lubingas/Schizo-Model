# app/equations.R
# ---------------------------
# The cost equations for the About tab, typeset with MathJax.
# The same equations, with a table of the inputs behind each symbol, are in the README;
# ip_costs(), op_costs(), residual_cost(), acute_cost() and third_line_cost() in R/functions.R follow them.
# ---------------------------

eq <- function(tex) div(class = "eq", paste0("\\[", tex, "\\]"))

cost_equations <- function() {
  withMathJax(
    p("Each year (cycle length \\(\\tau\\) = 1 year) a person's cost depends on the state they are in. Costs are discounted ",
      "and summed over the share of the cohort in each state in each cycle, \\(\\pi_{t,s}\\):"),
    eq(r"(\text{Cost} = \sum_{t=0}^{T-1} \frac{1}{(1+r)^{t}} \sum_{s} \pi_{t,s}\, C_s)"),
    h3("Side effects"),
    p("On drug \\(j\\), extrapyramidal symptoms (EPS) come from the odds ratio against placebo, and diabetes from the relative risk against haloperidol:"),
    eq(r"(\text{EPS}_j = \frac{\text{OR}_j\, p_0}{1 - p_0 + \text{OR}_j\, p_0} \qquad \text{DM}_j = p_{\text{DM}}\, \text{RR}_j)"),
    h3("Staff pay"),
    p("Hourly pay is annual salary over the hours worked in a year, \\(h_k = S_k / H\\), for psychiatrists (psy), psychiatric ",
      "clinical officers (co), psychiatric nurses (nu) and pharmacy assistants (ph)."),
    h3("A cycle in hospital"),
    p("The four relapse states: a stay of \\(L\\) days on chlorpromazine. Ward rounds are counted per week, so they are divided by 7 for a daily cost."),
    eq(r"(C^{\text{IP}}_{\text{staff}} = \left[\frac{1}{7}\left(R_{\text{maj}} \frac{d_{\text{maj}}}{60} \big(N^{\text{maj}}_{\text{psy}} h_{\text{psy}} + N^{\text{maj}}_{\text{co}} h_{\text{co}} + N^{\text{maj}}_{\text{nu}} h_{\text{nu}}\big) + R_{\text{min}} \frac{d_{\text{min}}}{60} \big(N^{\text{min}}_{\text{co}} h_{\text{co}} + N^{\text{min}}_{\text{nu}} h_{\text{nu}}\big)\right) + N_{\text{adm}} \frac{d_{\text{adm}}}{60} h_{\text{nu}}\right] \tau L)"),
    eq(r"(C^{\text{IP}}_{\text{lab}} = \frac{\text{CPI}_{2012}}{\text{CPI}_{2001}} \Big(c_{\text{FBC}} f_{\text{FBC}} + c_{\text{ESR}} f_{\text{ESR}} + p_{\text{TPHA}}\, c_{\text{TPHA}} f_{\text{TPHA}} + (1 - p_{\text{TPHA}})\, c_{\text{VDRL}} f_{\text{VDRL}} + c_{\text{HIV}} f_{\text{HIV}}\Big)\, 12\tau)"),
    eq(r"(C^{\text{IP}}_{\text{capital}} = k_{\text{IP}}\, \tau L \qquad C^{\text{IP}}_{\text{travel}} = 2 c_{\text{trip}}\, p_{\text{carer}} \frac{\tau L\, v_{\text{carer}}}{7} + 2 c_{\text{trip}} \qquad C^{\text{IP}}_{\text{work}} = \lambda (1 - u)\, \text{GDP}\, \tau)"),
    eq(r"(C^{\text{acute}} = L\, c_{\text{CPZ}} n_{\text{CPZ}} + \text{EPS}_{\text{CPZ}}\, L\, c_{\text{bz}} n_{\text{bz}} + \text{DM}_{\text{CPZ}}\, C_{\text{DM}} + C^{\text{IP}}_{\text{staff}} + C^{\text{IP}}_{\text{lab}} + C^{\text{IP}}_{\text{capital}} + C^{\text{IP}}_{\text{travel}} + C^{\text{IP}}_{\text{work}} + C_{\text{carer}})"),
    h3("A cycle as an outpatient"),
    p("\\(V\\) visits a month, so \\(12V\\tau\\) visits a cycle. The consultation is with a psychiatrist (a share \\(\\sigma_{\\text{psy}}\\) ",
      "of patients) or a clinical officer, and a nurse; then the drugs are dispensed."),
    eq(r"(C^{\text{OP}}_{\text{staff}} = \left[\frac{d_{\text{con}}}{60} \big(N_{\text{psy}}\, \sigma_{\text{psy}} h_{\text{psy}} + N_{\text{co}} (1 - \sigma_{\text{psy}}) h_{\text{co}} + N_{\text{nu}}\, \sigma_{\text{nu}} h_{\text{nu}}\big) + N_{\text{ph}} \frac{d_{\text{disp}}}{60} h_{\text{ph}}\right] 12 V \tau)"),
    eq(r"(C^{\text{OP}}_{\text{lab}} = \frac{\text{CPI}_{2012}}{\text{CPI}_{2001}} \big(c_{\text{FBC}} f_{\text{FBC}} + c_{\text{ESR}} f_{\text{ESR}}\big)\, 12\tau \qquad C^{\text{OP}}_{\text{capital}} = k_{\text{OP}}\, 12 V \tau)"),
    eq(r"(C^{\text{OP}}_{\text{travel}} = 2 c_{\text{trip}} (p_{\text{carer}} + 1)\, 12 V \tau \qquad C^{\text{OP}}_{\text{work}} = \lambda (1 - u) \frac{\text{GDP}}{264}\, 12 V \tau)"),
    p("On drug \\(j\\) (\\(n_j\\) tablets a day at \\(c_j\\) each; benzhexol for EPS for \\(\\delta_{\\text{EPS}}\\) years), and on fluphenazine decanoate (injections a month):"),
    eq(r"(C^{\text{res}}_j = 365\, c_j n_j + \text{EPS}_j\, 365\, c_{\text{bz}} n_{\text{bz}}\, \delta_{\text{EPS}} + \text{DM}_j\, C_{\text{DM}} + C^{\text{OP}}_{\text{staff}} + C^{\text{OP}}_{\text{lab}} + C^{\text{OP}}_{\text{capital}} + C^{\text{OP}}_{\text{travel}} + C^{\text{OP}}_{\text{work}} + C_{\text{carer}})"),
    eq(r"(C^{\text{res}}_{\text{FLU}} = 12\, c_{\text{FLU}} n_{\text{FLU}} + \text{EPS}_{\text{FLU}}\, 365\, c_{\text{bz}} n_{\text{bz}}\, \delta_{\text{EPS}} + p_{\text{DM,FLU}}\, C_{\text{DM}} + C^{\text{OP}}_{\text{staff}} + C^{\text{OP}}_{\text{lab}} + C^{\text{OP}}_{\text{capital}} + C^{\text{OP}}_{\text{travel}} + C^{\text{OP}}_{\text{work}} + C_{\text{carer}})"),
    h3("Carers"),
    p("Carers lose a share \\(\\ell\\) of their work time (264 work days a year), in hospital and as an outpatient. ",
      "Time off treatment (the stopped states) costs nothing."),
    eq(r"(C_{\text{carer}} = \lambda (1 - u)\, \text{GDP}\, \tau\, \ell\, \tfrac{264}{365})"),
    p(class = "small-print", "Every symbol is an input on the Inputs tab; the README on GitHub lists which input each one is.")
  )
}

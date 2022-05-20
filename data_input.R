library(readxl)

# load data
data = read_excel("data_input.xlsx", sheet = 1)
cpi = read_excel("data_input.xlsx", sheet = 2)
lifeTab = read_excel("data_input.xlsx", sheet = 3)
schizoTab = read_excel("data_input.xlsx", sheet = 4)
lifeExp = read_excel("data_input.xlsx", sheet = 5)

# colors
col <- brewer.pal(6,"Dark2")



# parameters for the one way sensitivity analysis


# inputs for PSA


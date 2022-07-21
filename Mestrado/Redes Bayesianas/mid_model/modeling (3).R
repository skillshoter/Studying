
setwd("D:/PhD/IA368 - Análise de dados Bayesiana_1_Sem_2022/Trabalho_1/codigo_R/mod7")

require(pacman)
p_load(tidyverse)
p_load(magrittr)
p_load(rstan)
p_load(ggplot2)
p_load(rstanarm)
p_load(glue)
p_load(scales)
p_load(gridExtra)

# read processed data file and remove world summarized information
df <- read.csv("processed_data_3.csv") %>%
  as_tibble() %>% 
  filter(iso_code != 'OWID_WRL')%>%
  mutate(intercept=1)

# target and deterministic column name
target_var <- "new_deaths"
x_given <- c("population", "group_idh")


# features column names
# Full model
#model_name <- "Full_model"
# x_names <- c(
#   "intercept",
#   "new_cases",
#   "hospital_beds_per_thousand",
#   "life_expectancy",
#   "median_age",
#   "diabetes_prevalence",
#   "cardiovasc_death_rate",
#   "population_density",
#   "gdp_per_capita"
# )
# Reduced model
model_name <- "Reduced_model"
x_names <- c(
  "intercept",
  "new_cases",
  "life_expectancy"
)


# column names to keep
columns_to_keep <- c(target_var, x_given, x_names)
df <- df %>% select(all_of(columns_to_keep))

# number of hierarchical levels
n_idh_groups <- df %>%
  pull(group_idh) %>%
  unique() %>% 
  length()

# force population column to be integer
df$population <- df %>%
  mutate(population = as.integer(round(population, 0))) %>%
  pull(population)

# create data for stan model
covid_data <- list(
  N = nrow(df),                             # quantidade de paises
  D = length(x_names),                      # número de preditores
  L = n_idh_groups,                         # quantidade de levels da hierarquia
  y = df %>% pull(new_deaths),              # número de mortes em cada pais
  pop = df %>% pull(population),            # população em cada pais
  ll = df %>% pull(group_idh),              # identificador da classe da hierarquia
  x = as.matrix(                            # covariáveis 
    df %>%
      select(all_of(x_names))
    )
)

# run HMC
n_chain <- 4
n_iter <- 2000
fit <- stan(
  file = 'stan_code.stan',
  data = covid_data,
  chains = n_chain,
  control = list(max_treedepth = 12),
  iter = n_iter
)

# define all parameters present in the model
param_names <- c("beta", "mu", "sigma", "y_rep")

# extract all simulated values for the param_names, after the initial warmup
sim_chains <- rstan::extract(
  fit,
  pars = param_names,
  permuted = TRUE,
  inc_warmup = FALSE,
  include = TRUE
)

curr_param_group <- "beta"
plot_type <- "hist"

# simulations of the current group of parameter
group_sim <- sim_chains[[curr_param_group]]

tb_sim <- group_sim %>% as_tibble() %>% gather() %>%
  mutate(
    key = gsub("V", "1.", key)
  ) %>% 
  separate(key, c("hierarchical_group", "feature"), remove = FALSE) %>% 
  group_by(feature, hierarchical_group) %>% 
  mutate(n_count = 1:n()) %>% 
  ungroup() %>% 
  mutate(chain = cut(n_count, breaks = seq(0, n_iter*n_chain/2, n_iter/2), labels = 1:n_chain)) %>% 
  group_by(feature, hierarchical_group, chain) %>% 
  mutate(n_count = 1:n()) %>% 
  ungroup() %>% 
  mutate(
    feature_name = case_when(
      feature == 1 ~ x_names[1],
      feature == 2 ~ x_names[2],
      feature == 3 ~ x_names[3],
      feature == 4 ~ x_names[4],
      feature == 5 ~ x_names[5],
      feature == 6 ~ x_names[6],
      TRUE ~ "not_mapped"
    ),
    idh_name = paste("IDH:", hierarchical_group)
  )

unique_plot <- tb_sim %>% select(feature_name, idh_name) %>% unique()
curr_n_h <- tb_sim %>% pull(hierarchical_group) %>% unique() %>% length()
curr_n_f <- tb_sim %>% pull(feature_name) %>% unique() %>% length()
p <- list()
for (i in 1:nrow(unique_plot)) {
  curr_row <- unique_plot[i, ]
  curr_tb <- tb_sim %>%
    filter(feature_name == curr_row$feature_name & idh_name == curr_row$idh_name)
  
  if (plot_type == "chains") {
    curr_plot <- ggplot(curr_tb, aes(x = n_count, y = value, color = chain)) +
      geom_line()
    x_lab = "Iteration"
    y_lab = "Value"
  } else {
    curr_plot <- ggplot(curr_tb, aes(x = value, color = chain)) +
      geom_density()
    x_lab = "Value"
    y_lab = "Density"
  }
  if (curr_n_h == 1) {
    title_ <- glue("{curr_param_group} | {curr_row$feature_name}")
  } else {
    title_ <- glue("{curr_param_group} | {curr_row$feature_name} | {curr_row$idh_name}")
  }
  curr_plot <- curr_plot + labs(title = title_, x = x_lab, y = y_lab) +
    theme_bw() + 
    theme(
      text = element_text(size=8),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
  p[[i]] <- curr_plot
}

pdf(
  file = glue("{curr_param_group}_{plot_type}.pdf"),
  width = 15,
  height = 30
)
do.call("grid.arrange", c(p, ncol = curr_n_h, nrow = curr_n_f))
dev.off()

print(fit)
pdf(
  file = glue("params_summary.pdf"),
  width = 10,
  height = 10
)
plot(fit, pars = c("mu", "sigma", "beta"))
dev.off()

# pairs(fit, pars = c("beta"))
# matrix_of_draws <- as.matrix(fit)

write.csv(sim_chains,glue("{model_name}_chains.csv"),row.names = FALSE)
write.csv(summary(fit)$summary,glue('{model_name}_summary.csv'))
write.csv(x_names,glue("{model_name}_x_names.csv"))
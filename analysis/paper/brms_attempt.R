library(tidyverse)
library(tidybayes)
library(here)
library(brms)

peri_all <- read_csv(here("analysis/data/derived_data/FA_peri_sums_areal.csv"))


FA_groups_areal <- peri_all|>
  rename_with(~str_remove(.,"_mean"))|>
  rowwise()|>
  mutate(Tot_FA = sum(across(all_of(4:40))),
         MUFA = sum(c_across(contains(":1"))),
         MUFA_16 = sum(c_across(contains("16:1"))),
         MUFA_18 = sum(c_across(contains("18:1"))),
         SAFA = sum(c_across(c(matches(":0"),-matches("i-")))),
         LSAFA = sum(c_across(c(`22:0`, `24:0`,`23:0`))),
         brFA = sum(c_across(contains("i-"))),
         diatom_16PUFA = sum(c_across(c(`16:2w7`,`16:3w4`))),
         green_16PUFA = sum(c_across(c(`16:3w3`))),
         w3_sum = sum(c_across(contains("w3"))),
         w6_sum = sum(c_across(contains("w6"))),
         w3_w6 = w3_sum/w6_sum, #note that the ratio of two random variables does not have a well-defined variance/sd, or at least not a closed form solution
         EPA = `20:5w3`,
         DHA = `22:6w3`,
         DPA = `22:5w3`,
         LIN = `18:2w6`,
         SDA = `18:4w3`,
         ALA = `18:3w3`,
         ARA = `20:4w6`,
         GLA = `18:3w6`,
         sum_PUFA = sum(across(c(EPA,DHA,DPA,LIN,SDA,LIN,ARA,GLA))),
         #sum variances and sqrt for sd
         Tot_FA_sd = sqrt(sum(across(contains("_var")))),
         MUFA_sd = sqrt(sum(c_across(contains(":1")&contains("_var")))),
         MUFA_16_sd = sqrt(sum(c_across(contains("16:1")&contains("_var")))),
         MUFA_18_sd = sqrt(sum(c_across(contains("18:1")&contains("_var")))),
         SAFA_sd = sqrt(sum(c_across(c(matches(":0"),-matches("i-"))))),
         LSAFA_sd = sqrt(sum(c_across(c(`22:0_var`, `24:0_var`,`23:0_var`)))),
         brFA_sd = sqrt(sum(c_across(contains("i-")&contains("_var")))),
         diatom_16PUFA_sd = sqrt(sum(c_across(c(`16:2w7_var`,`16:3w4_var`)))),
         green_16PUFA_sd = sqrt(sum(c_across(contains("16:3w3")&contains("_var")))),
         w3_sum = sum(c_across(contains("w3")&contains("_var"))),
         w6_sum = sum(c_across(contains("w6")&contains("_var"))),
         EPA_sd = sqrt(`20:5w3_var`),
         DHA_sd = sqrt(`22:6w3_var`),
         DPA_sd = sqrt(`22:5w3_var`),
         LIN_sd = sqrt(`18:2w6_var`),
         SDA_sd = sqrt(`18:4w3_var`),
         ALA_sd = sqrt(`18:3w3_var`),
         ARA_sd = sqrt(`20:4w6_var`),
         GLA_sd = sqrt(`18:3w6_var`),
         sum_PUFA = sum(across(c(EPA,DHA,DPA,LIN,SDA,LIN,ARA,GLA))),
         .keep = "unused", .after = 1)



#brms model for EPA
d <- FA_groups_areal|>
  ungroup()|>
  select(Site, LSAFA, LSAFA_sd, Canopy_Cover_per, NO3_ug.L, EcoRegion)|>
  mutate(LSAFA_obs = (LSAFA-mean(LSAFA))/LSAFA_sd,
         LSAFA_sd_std = (LSAFA)/sd(LSAFA),
         CC_std = Canopy_Cover_per - mean(Canopy_Cover_per)/sd(Canopy_Cover_per),
         NO3_std = NO3_ug.L - mean(NO3_ug.L)/sd(NO3_ug.L))

dlist <- list(
  LSAFA_obs = d$LSAFA_obs,
  LSAFA_sd  = d$LSAFA_sd_std,
  CC     = d$CC_std,
  NO3     = d$NO3_std,
  EcoRegion = d$EcoRegion)


b15.1 <-
  brm(data = dlist,
      family = gaussian,
      EPA_obs | mi(EPA_sd) ~ 1 + CC + NO3,
      prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(exponential(1), class = sigma)),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 15,
      save_pars = save_pars(latent = TRUE)) #save_pars super important to getting everything to pop up

posterior_summary(b15.1)

#no error incorporated
b15.2 <-
  brm(data = dlist,
      family = gaussian,
      EPA_obs  ~ 1 + CC + NO3,
      prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(exponential(1), class = sigma)),
      iter = 2000, warmup = 1000, cores = 4, chains = 4,
      seed = 15,
      save_pars = save_pars(latent = TRUE)) #save_pars super important to getting everything to pop up

posterior_summary(b15.2)

b15.3 <-
  brm(data = dlist,
      family = gaussian,
      LSAFA_obs | mi(LSAFA_sd) ~ 1 + CC+NO3+EcoRegion,
      prior = c(prior(normal(0, 1), class = Intercept),
                prior(normal(0, 1), class = b),
                prior(exponential(1), class = sigma)),
      iter = 4000, warmup = 1000, cores = 4, chains = 4,
      seed = 15,
      save_pars = save_pars(latent = TRUE)) #save_pars super important to getting everything to pop up

posterior_summary(b15.3)



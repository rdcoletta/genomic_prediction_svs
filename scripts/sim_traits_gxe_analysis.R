library(data.table)
library(ggplot2)
library(dplyr)

# install.packages('reshape', repos='http://cran.rstudio.com/')
# install.packages("plantbreeding", repos="http://R-Forge.R-project.org")
# 
# library(plantbreeding)

# install.packages("metan")

library(metan)


qtn_number <- 100
h2 <- 0.5
effect <- 0.1
# mean_cor <- 0.35
mean_cor <- 0.8

sim_folder <- paste0("analysis/trait_sim_mult-env/additive_model/test/genetic_cor_matrix/", 
                     qtn_number, "-QTNs_from_SNP/",
                     h2, "-heritability/", 
                     "effect", effect, "/",
                     "mean-gen-cor_", mean_cor)

# load sim traits
sim_traits <- paste0(sim_folder, "/simulated_traits_summary.txt")
sim_traits <- fread(sim_traits, header = TRUE, data.table = FALSE)

# get trait means for each genotype in each environment
sim_traits_mean <- sim_traits %>% 
  dplyr::group_by(genotype, environment) %>% 
  dplyr::summarize(trait_mean = mean(trait_value))

# get environmental mean
sim_envs_mean <- sim_traits %>% 
  dplyr::group_by(environment) %>% 
  dplyr::summarize(env_mean = mean(trait_value)) %>% 
  dplyr::arrange(env_mean)

# reorder factor based on environmental mean
sim_traits_mean$environment <- factor(sim_traits_mean$environment, 
                                      levels = sim_envs_mean$environment)


ggplot(sim_traits_mean) +
  geom_line(aes(x = environment, y = trait_mean, group = genotype), alpha = 0.2) +
  geom_line(data = sim_envs_mean, aes(x = environment, y = env_mean, group = 1), color = "red")

ggsave(filename = paste0(sim_folder, "/mean_genotypes_by_envs.pdf"), device = "pdf")

# sim_traits %>% 
#   filter(genotype == "B73*LH82-B-B-10-1-1-B-B") %>% 
#   ggplot(aes(x = environment, y = trait_value, group = genotype)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE)
# 
# 
# 
# sim_traits_B73xLH82 <- sim_traits[grep("B73*LH82", sim_traits$genotype, fixed = TRUE), ]
# sim_traits_B73xLH82$genotype <- factor(sim_traits_B73xLH82$genotype)
# sim_traits_B73xLH82$environment <- factor(sim_traits_B73xLH82$environment)
# sim_traits_B73xLH82$rep <- factor(sim_traits_B73xLH82$rep)
# 
# 
# 
# reg <- sim_traits_B73xLH82 %>% 
#   ge_reg(environment, genotype, rep,
#          resp = trait_value)
# print(reg)
# 
# reg$trait_value$data %>%
#   ggplot(aes(x = IndAmb, y = Y, color = GEN)) +
#   geom_point(show.legend = FALSE) +
#   geom_smooth(method = "lm", se = FALSE, show.legend = FALSE)
# 
# reg$trait_value$data %>% 
#   group_by(GEN) %>% 
#   do({
#     mod = lm(Y ~ IndAmb, data = .)
#     data.frame(Intercept = coef(mod)[1],
#                Slope = coef(mod)[2])
#   })
# 
# 
# 
# reg$trait_value$data %>%
#   filter(IndAmb < 0.04) %>% 
#   ggplot(aes(x = IndAmb, y = Y, color = GEN)) +
#   geom_point(show.legend = FALSE) +
#   geom_smooth(method = "lm", se = FALSE, show.legend = FALSE)
# 
# reg$trait_value$data %>% 
#   filter(IndAmb < 0.04) %>% 
#   group_by(GEN) %>% 
#   do({
#     mod = lm(Y ~ IndAmb, data = .)
#     data.frame(Intercept = coef(mod)[1],
#                Slope = coef(mod)[2])
#   })






sim_traits$genotype <- factor(sim_traits$genotype)
sim_traits$environment <- factor(sim_traits$environment)
sim_traits$rep <- factor(sim_traits$rep)

reg_full <- sim_traits %>% 
  ge_reg(environment, genotype, rep,
         resp = trait_value)

View(reg_full$trait_value$data)
View(reg_full$trait_value$anova)
View(reg_full$trait_value$regression)

reg_full$trait_value$data %>%
  ggplot(aes(x = IndAmb, y = Y, group = GEN)) +
  # geom_point(show.legend = FALSE) +
  geom_line(stat = "smooth", method = "lm", formula = y ~ x, color = "black", alpha = 0.5, se = FALSE)




high_slope <- reg_full$trait_value$data %>% 
  group_by(GEN) %>% 
  do({
    mod = lm(Y ~ IndAmb, data = .)
    data.frame(Intercept = coef(mod)[1],
               Slope = coef(mod)[2])
  }) %>% 
  filter(Slope > 70) %>% 
  pull(GEN)

high_slope <- as.character(high_slope)

reg_full$trait_value$data %>%
  filter(GEN %in% high_slope) %>% 
  ggplot(aes(x = IndAmb, y = Y, color = GEN)) +
  geom_point(show.legend = FALSE) +
  geom_smooth(method = "lm", se = FALSE, show.legend = FALSE)



high_R2 <- reg_full$trait_value$regression %>% 
  filter(R2 > 0.5)



# high slope, but also good fit
high_slope <- reg_full$trait_value$regression %>% 
  filter(slope > 1.5 | slope < -1.5) %>% 
  filter(R2 > 0.5) %>% 
  pull(GEN)

reg_full$trait_value$data$family <- sapply(as.character(reg_full$trait_value$data$GEN),
                                           function(x) unlist(strsplit(x, split = "-", fixed = TRUE))[1])

reg_full$trait_value$data %>% 
ggplot(aes(x = IndAmb, y = Y, group = GEN)) +
  # geom_point(show.legend = FALSE) +
  geom_line(stat = "smooth", method = "lm", formula = y ~ x, color = "black", alpha = 0.5, se = FALSE) +
  geom_line(data = reg_full$trait_value$data[which(reg_full$trait_value$data$GEN %in% high_slope), ],
            aes(x = IndAmb, y = Y, group = GEN),
            stat = "smooth", method = "lm", formula = y ~ x,
            color = "red", alpha = 1, se = FALSE) +
  facet_wrap(~family, nrow = 2) +
  labs(x = "environment index",
       y = "simulated trait value") +
  theme(axis.text.x = element_text(size = 8))
  
ggsave(filename = paste0(sim_folder, "/finlay-wilkinson_stability_per_family.pdf"), device = "pdf",
       width = 14, height = 7)




install.packages("statgenGxE")
library(statgenGxE)

simTD <- statgenSTA::createTD(data = sim_traits, genotype = "genotype", trial = "environment")

plot(simTD, plotType = "box", traits = "trait_value", orderBy = "descending", colorTrialBy = "trial", colorGenoBy = "genotype")
plot(simTD, plotType = "scatter", traits = "trait_value")

simFW <- gxeFw(TD = simTD, trait = "trait_value", maxIter = 100)
summary(simFW)

plot(simFW, plotType = "scatter")
# plot(simFW, plotType = "scatter", colorGenoBy = "geneticGroup")
# need to create a column in sim_traits with "geneticGroup" as the family name of a RIL

plot(simFW, plotType = "line")
# plot(dropsFW, plotType = "line", colorGenoBy = "geneticGroup")

# https://cran.r-project.org/web/packages/statgenGxE/vignettes/statgenGxE.html
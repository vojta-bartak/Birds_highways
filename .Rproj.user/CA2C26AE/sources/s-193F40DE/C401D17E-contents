library(tidyverse)
library(readxl)
library(sf)
library(ncf)
library(lme4)
library(lmerTest)
library(performance)
library(sjPlot)
library(spdep)

# Data -------------------------------------------------------------------------------------------------------------------------
df.old <- read_xlsx("Data_DP.xlsx", sheet=2) %>% 
  #mutate_at(names(.)[is.character(names(.)) & !grepl("souradnice",names(.))], as.factor) %>% 
  mutate(lat = sapply(souradnice, function(sour) {strsplit(sour, split="N, ")[[1]][1] %>% as.numeric}),
         lon = sapply(souradnice, function(sour) {
           strsplit(sour, split="N, ")[[1]][2] %>% strsplit("E") %>% "[["(1) %>% "["(1) %>% as.numeric}),
         transekt = gsub('[[:alpha:]]', '', KOD_Plochy),
         forest_sc = scale(forest),
         village_sc = scale(village),
         water_sc = scale(water),
         distance_sc = scale(distance),
         ruderal = as.factor(ruderal),
         louka= as.factor(louka),
         field = as.factor(field)
         )
df.sf <- df.old %>% 
  st_as_sf(coords=c("lon","lat"), crs = 4326) %>% 
  st_transform(5514)
plot(df.sf %>% st_geometry)
plot(df.sf %>% filter(startsWith(KOD_Plochy, "12")) %>% st_geometry)
birds <- read_xlsx("Data_DP.xlsx", sheet=3)

library(tmap)
tmap_mode("view")
tm_shape(df.sf %>% filter(dalnice=="D10")) + tm_dots("transekt")

# Spatial autocorrelation ------------------------------------------------------------------------------------------------------
df.cor <- df.sf %>% 
  st_coordinates %>% 
  cbind(df$Pocet_druhu)

correlog(df.cor[,1], df.cor[,2], df.cor[,3], increment = 500, resamp = 100) %>% plot(xlim = c(0,10000), ylim=c(-1,1))

# Models -----------------------------------------------------------------------------------------------------------------------

# full model
m <- glmer(Pocet_druhu~distance_sc+E2+E3+forest_sc+village_sc+water_sc+
             (1|pozorovatel/transekt), data=df, family="poisson")

# random effects selection
anova(m, update(m, ~.-(1|pozorovatel/transekt)+(1|pozorovatel)))
anova(update(m, ~.-(1|pozorovatel/transekt)+(1|pozorovatel)), 
      glm(Pocet_druhu~E2+E3+forest_sc+village_sc+water_sc, data=df, family="poisson"))
m <- glmer(Pocet_druhu~E2+E3+forest_sc+village_sc+water_sc+distance_sc+
             (1|pozorovatel/transekt), data=df, family="poisson")
car::vif(m)

# fixed effects selection
summary(m)
drop1(m, test="Chisq")
drop1(update(m, ~.-forest_sc))
drop1(update(m, ~.-forest_sc-E2))
m <- glmer(Pocet_druhu~E2+E3+village_sc+water_sc+(1|pozorovatel), data=df, family="poisson")

# overdispersion
check_overdispersion(m)

# model outputs
summary(m)
plot_model(m, type="eff", terms = "E3", show.data = T)
eff_plots <- plot_model(m, type="eff", show.data = T)
eff_plots$E3
eff_plots$village_sc
eff_plots$water_sc
eff_plots$distance_sc
plot_model(m, type="re")
(rsq <- r2(m))
(importances <- sapply(names(eff_plots), function(pr) rsq$R2_marginal - r2(update(m, paste("~.-", pr, sep="")))$R2_marginal) %>% 
  set_names(names(eff_plots)))
barplot(importances %>% sort(decreasing=T))


# farmland, Woodland, hnízdní gildy
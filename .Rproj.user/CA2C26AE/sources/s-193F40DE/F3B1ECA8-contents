library(tidyverse)
library(readxl)
library(lme4)

# Data -------------------------------------------------------------------------------------------------------------------------
df <- read_xlsx("data.xlsx", sheet=1) %>% 
  mutate(forest_sc = scale(Forest),
         village_sc = scale(Village),
         distance_sc = scale(Distance),
         traffic_sc = scale(Traffic_volume),
         alt_sc = scale(Altitude),
         e2_sc = scale(E2),
         e3_sc = scale(E3),
         water = as.factor(Water),
         ruderal = as.factor(Ruderal),
         meadow = as.factor(Meadow),
         field = as.factor(Field),
         noise = as.factor(Noise),
         transect = sapply(ID, function(id) substr(id, 1, str_length(id)-1))
  )
orig.preds <- c(forest_sc="Forest", village_sc="Village", distance_sc="Distance", traffic_sc="Traffic_volume", alt_sc="Altitude",
                e2_sc="E2", e3_sc="E3", water="Water", ruderal="Ruderal", meadow="Meadow", field="Field", noise="Noise")

stepAIC <- function(model, thres = 2){
  models <- list(model)
  dr <- drop1(model)
  aics <- c(dr[1,]$AIC)
  print(paste("AIC = ", aics[length(aics)]))
  dr_var <- rownames(dr[dr$AIC==min(dr$AIC[-1]),])
  dr_vars <- c("<none>", dr_var)
  while (nrow(dr) > 1) {
    print(paste("Droping", dr_var))
    model <- update(model, paste("~.-", dr_var, sep=""))
    dr <- drop1(model)
    aics <- c(aics, dr[1,]$AIC)
    print(paste("AIC = ", aics[length(aics)]))
    dr_var <- rownames(dr[dr$AIC==min(dr$AIC[-1]),])
    dr_vars <- c(dr_vars, dr_var)
    models <- c(models, list(model))
  }
  plot(aics~c(1:length(dr_vars)), type="l", xlab="", ylab="AIC", xaxt="n")
  points(aics~c(1:length(dr_vars)))
  points(min(aics)~which(aics==min(aics)), col = "red", pch = 16)
  axis(1, at = 1:length(aics), labels = dr_vars, las=3)
  return(models[[which(aics == min(aics))]])
}

# Modeling strategies ----------------------------------------------------
m.Species1 <- glmer(Species~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                      (1|Highway), data=df, family="poisson")
m.Species2 <- glmer(Species~noise*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+alt_sc)+
                      (1|Highway), data=df, family="poisson")
m.Species3 <- glmer(Species~noise+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+alt_sc+
                      (1|Highway), data=df, family="poisson")
AIC(m.Species1, m.Species2, m.Species3)

m.Species1.sel <- stepAIC(m.Species1)
m.Species2.sel <- stepAIC(m.Species2)
m.Species3.sel <- stepAIC(m.Species3)

AIC(m.Species1.sel, m.Species2.sel, m.Species3.sel)

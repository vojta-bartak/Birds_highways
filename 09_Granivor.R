library(tidyverse)
library(readxl)
library(sf)
library(ncf)
library(lme4)
library(lmerTest)
library(performance)
library(sjPlot)
library(car)

# Data -------------------------------------------------------------------------------------------------------------------------
df <- read_xlsx("Data_Štěpán_Hladík 30.8.2023.xlsx", sheet=1) %>% 
  mutate(forest_sc = scale(Forest) %>% c,
         village_sc = scale(Village) %>% c,
         distance_sc = scale(Distance) %>% c,
         traffic_sc = scale(Traffic_volume) %>% c,
         alt_sc = scale(Altitude) %>% c,
         e2_sc = scale(E2) %>% c,
         e3_sc = scale(E3) %>% c,
         water = as.factor(Water),
         ruderal = as.factor(Ruderal),
         meadow = as.factor(Meadow),
         field = as.factor(Field),
         transect = sapply(ID, function(id) substr(id, 1, str_length(id)-1))
  )
orig.preds <- c(forest_sc="Forest", village_sc="Village", distance_sc="Distance", traffic_sc="Traffic_volume", alt_sc="Altitude",
                e2_sc="E2", e3_sc="E3", water="Water", ruderal="Ruderal", meadow="Meadow", field="Field")

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


# Granivor richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Granivor <- glmer(Granivor~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                   (1|Highway/transect), data=df, family="poisson")
vif(m.Granivor)

# random effects selection
m2.Granivor <- update(m.Granivor, ~.-(1|Highway/transect)+(1|Highway))
m3.Granivor <- glm(Granivor~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc, data=df, family="poisson")
anova(m.Granivor, m2.Granivor, m3.Granivor)
m.Granivor <- glm(Granivor~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc, data=df, family="poisson")

# fixed effects selection
summary(m.Granivor)
m.Granivor <- stepAIC(m.Granivor)

# diagnostics
check_overdispersion(m.Granivor)

# model outputs -------------------
(s.Granivor <- summary(m.Granivor))
(drop.Granivor <- drop1(m.Granivor, test="Chisq"))
(r2.Granivor <- r2(m.Granivor))

(preds.Granivor <- all.vars(formula(m.Granivor))[-1][-10])
preds.cont.Granivor <- preds.Granivor[sapply(preds.Granivor, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Granivor <- preds.Granivor[sapply(preds.Granivor, function(pred) is.factor(pull(df, pred)))]
# varimp.r2.Granivor <- sapply(preds.Granivor, function(pred){
#   r2.Granivor$R2_marginal - r2(update(m.Granivor, paste("~.-", pred, sep="")))$R2_marginal
# }) %>% set_names(preds.Granivor)
varimp.r2.Granivor <- sapply(preds.Granivor, function(pred){
  r2.Granivor$R2_Nagelkerke - r2(update(m.Granivor, paste("~.-", pred, sep="")))$R2_Nagelkerke
}) %>% set_names(preds.Granivor)

barplot(sort(varimp.r2.Granivor, decreasing = T), las=3)

# predictions
nds.cont.Granivor <- lapply(preds.cont.Granivor, function(pred){
  nd <- lapply(preds.cont.Granivor, function(pr) rep(mean(pull(df, pr)), 100)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Granivor)
  nd[,pred] <- seq(min(pull(df, pred)),  max(pull(df, pred)), l=100)
  mm <- model.matrix(as.formula(paste("~", paste(preds.Granivor, collapse = "+"), sep="")), data=nd)
  nd$y <- predict(m.Granivor, newdata = nd, type = "link", se=F)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Granivor), mm))
  nd
}) %>% set_names(preds.cont.Granivor)

# prediction plots
plots.cont.Granivor <- lapply(preds.cont.Granivor, function(pred){
  nd <- nds.cont.Granivor[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Granivor), alpha = .1) +
    geom_line() +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar))), alpha= .3) +
    labs(y="Granivor richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Granivor)

cowplot::plot_grid(plotlist=plots.cont.Granivor)
ggsave("prediction.plots.Granivor.png", 
       plot = cowplot::plot_grid(plotlist=plots.cont.Granivor), 
       width = 27, height = 18, units = "cm", dpi=300)

outs.Granivor = list(
  model = m.Granivor,
  summary = s.Granivor,
  drop = drop.Granivor,
  r2 = r2.Granivor,
  varimp = varimp.r2.Granivor,
  preds = preds.Granivor,
  preds.cont = preds.cont.Granivor,
  preds.fact = NULL,
  nds.cont = nds.cont.Granivor,
  nds.fact = NULL
)
save(outs.Granivor, file = "outs.Granivor.RData")

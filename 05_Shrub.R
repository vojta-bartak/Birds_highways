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


# Shrub richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Shrub <- glmer(Shrub~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                    (1|Highway/transect), data=df, family="poisson")
vif(m.Shrub)

# random effects selection
m2.Shrub <- update(m.Shrub, ~.-(1|Highway/transect)+(1|Highway))
m3.Shrub <- glm(Shrub~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc, data=df, family="poisson")
anova(m.Shrub, m2.Shrub, m3.Shrub)
m.Shrub <- glm(Shrub~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc, data=df, family="poisson")

# fixed effects selection
summary(m.Shrub)
m.Shrub <- stepAIC(m.Shrub)

# diagnostics
check_overdispersion(m.Shrub)

# model outputs -------------------
(s.Shrub <- summary(m.Shrub))
(drop.Shrub <- drop1(m.Shrub, test="Chisq"))
(r2.Shrub <- r2(m.Shrub))

(preds.Shrub <- all.vars(formula(m.Shrub))[-1][-10])
preds.cont.Shrub <- preds.Shrub[sapply(preds.Shrub, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Shrub <- preds.Shrub[sapply(preds.Shrub, function(pred) is.factor(pull(df, pred)))]
# varimp.r2.Shrub <- sapply(preds.Shrub, function(pred){
#   r2.Shrub$R2_marginal - r2(update(m.Shrub, paste("~.-", pred, sep="")))$R2_marginal
# }) %>% set_names(preds.Shrub)
varimp.r2.Shrub <- sapply(preds.Shrub, function(pred){
  r2.Shrub$R2_Nagelkerke - r2(update(m.Shrub, paste("~.-", pred, sep="")))$R2_Nagelkerke
}) %>% set_names(preds.Shrub)

barplot(sort(varimp.r2.Shrub, decreasing = T), las=3)

# predictions
nds.cont.Shrub <- lapply(preds.cont.Shrub, function(pred){
  nd <- lapply(preds.cont.Shrub, function(pr) rep(mean(pull(df, pr)), 100)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Shrub) %>% 
    cbind(lapply(preds.fact.Shrub, function(pr) factor(rep(levels(pull(df, pr))[1], 100), levels=levels(pull(df, pr)))) %>% 
            bind_cols %>% 
            set_names(preds.fact.Shrub))
  nd[,pred] <- seq(min(pull(df, pred)),  max(pull(df, pred)), l=100)
  mm <- model.matrix(as.formula(paste("~", paste(preds.Shrub, collapse = "+"), sep="")), data=nd)
  nd$y <- predict(m.Shrub, newdata = nd, type = "link", se=F)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Shrub), mm))
  nd
}) %>% set_names(preds.cont.Shrub)

nds.fact.Shrub <- lapply(preds.fact.Shrub, function(pred){
  nd <- lapply(preds.cont.Shrub, function(pr) rep(mean(pull(df, pr)), length(levels(pull(df, pred))))) %>% 
    bind_cols %>% 
    set_names(preds.cont.Shrub) %>% 
    cbind(lapply(preds.fact.Shrub, function(pr) 
      factor(rep(levels(pull(df, pr))[1], length(levels(pull(df, pred)))), levels=levels(pull(df, pr)))) %>% 
        bind_cols %>% 
        set_names(preds.fact.Shrub))
  nd[,pred] <- factor(levels(pull(df, pred)), levels = levels(pull(df, pred)))
  mm <- model.matrix(as.formula(paste("~", paste(preds.Shrub, collapse = "+"), sep="")), data=nd)
  nd$y <- predict(m.Shrub, newdata = nd, type = "link", se=F)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Shrub), mm))
  nd
}) %>% set_names(preds.fact.Shrub)

# prediction plots
plots.cont.Shrub <- lapply(preds.cont.Shrub, function(pred){
  nd <- nds.cont.Shrub[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Shrub), alpha = .1) +
    geom_line() +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar))), alpha= .3) +
    labs(y="Shrub richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Shrub)
plots.fact.Shrub <- lapply(preds.fact.Shrub, function(pred){
  nd <- nds.fact.Shrub[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred))) +
    geom_jitter(data=df, aes(y=Shrub), alpha = .1, width=.2) +
    geom_point() +
    geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                      ymax = exp(y + 1.96*sqrt(pvar))),
                  width=.2) +
    labs(y="Shrub richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.fact.Shrub)

cowplot::plot_grid(plotlist=c(plots.cont.Shrub, plots.fact.Shrub))
ggsave("prediction.plots.Shrub.png", 
       plot = cowplot::plot_grid(plotlist=c(plots.cont.Shrub, plots.fact.Shrub)), 
       width = 27, height = 18, units = "cm", dpi=300)

outs.Shrub = list(
  model = m.Shrub,
  summary = s.Shrub,
  drop = drop.Shrub,
  r2 = r2.Shrub,
  varimp = varimp.r2.Shrub,
  preds = preds.Shrub,
  preds.cont = preds.cont.Shrub,
  preds.fact = preds.fact.Shrub,
  nds.cont = nds.cont.Shrub,
  nds.fact = nds.fact.Shrub
)
save(outs.Shrub, file = "outs.Shrub.RData")

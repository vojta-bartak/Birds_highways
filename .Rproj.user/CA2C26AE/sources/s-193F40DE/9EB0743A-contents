library(tidyverse)
library(readxl)
library(sf)
library(ncf)
library(lme4)
library(lmerTest)
library(performance)
library(sjPlot)
library(spdep)
library(car)

# Data -------------------------------------------------------------------------------------------------------------------------
df <- read_xlsx("Data_Štěpán_Hladík 30.8.2023.xlsx", sheet=1) %>% 
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


# Species richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Species <- glmer(Species~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                     (1|Highway/transect), data=df, family="poisson")
vif(m)

# random effects selection
m2.Species <- update(m.Species, ~.-(1|Highway/transect)+(1|Highway))
m3.Species <- glm(Species~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc, data=df, family="poisson")
anova(m.Species, m2.Species, m3.Species)
m.Species <- glmer(Species~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                     (1|Highway), data=df, family="poisson")

# fixed effects selection
summary(m.Species)
m.Species <- stepAIC(m.Species)

# diagnostics
check_overdispersion(m.Species)

# model outputs -------------------
(s.Species <- summary(m))
(drop.Species <- drop1(m.Species, test="Chisq"))
(r2.Species <- r2(m.Species))

preds.Species <- all.vars(formula(m.Species))[-1][-10]
preds.cont.Species <- preds.Species[sapply(preds.Species, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Species <- preds.Species[sapply(preds.Species, function(pred) is.factor(pull(df, pred)))]
varimp.r2.Species <- sapply(preds.Species, function(pred){
  r2.Species$R2_marginal - r2(update(m.Species, paste("~.-", pred, sep="")))$R2_marginal
}) %>% set_names(preds.Species)
barplot(sort(varimp.r2.Species, decreasing = T), las=3)

# predictions
nds.cont.Species <- lapply(preds.cont.Species, function(pred){
  nd <- lapply(preds.cont.Species, function(pr) rep(mean(pull(df, pr)), 100)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Species) %>% 
    cbind(lapply(preds.fact.Species, function(pr) factor(rep(levels(pull(df, pr))[1], 100), levels=levels(pull(df, pr)))) %>% 
            bind_cols %>% 
            set_names(preds.fact.Species))
  nd[,pred] <- seq(min(pull(df, pred)),  max(pull(df, pred)), l=100)
  mm <- model.matrix(as.formula(paste("~", paste(preds.Species, collapse = "+"), sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Species)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Species), mm))
  nd
}) %>% set_names(preds.cont.Species)

nds.fact.Species <- lapply(preds.fact.Species, function(pred){
  nd <- lapply(preds.cont.Species, function(pr) rep(mean(pull(df, pr)), length(levels(pull(df, pred))))) %>% 
    bind_cols %>% 
    set_names(preds.cont.Species) %>% 
    cbind(lapply(preds.fact.Species, function(pr) 
      factor(rep(levels(pull(df, pr))[1], length(levels(pull(df, pred)))), levels=levels(pull(df, pr)))) %>% 
            bind_cols %>% 
            set_names(preds.fact.Species))
  nd[,pred] <- factor(levels(pull(df, pred)), levels = levels(pull(df, pred)))
  mm <- model.matrix(as.formula(paste("~", paste(preds.Species, collapse = "+"), sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Species)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Species), mm))
  nd
}) %>% set_names(preds.fact.Species)

# prediction plots
plots.cont.Species <- lapply(preds.cont.Species, function(pred){
  nd <- nds.cont.Species[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Species), alpha = .1) +
    geom_line() +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar))), alpha= .3) +
    labs(y="Species richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Species)
plots.fact.Species <- lapply(preds.fact.Species, function(pred){
  nd <- nds.fact.Species[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred))) +
    geom_jitter(data=df, aes(y=Species), alpha = .1, width=.2) +
    geom_point() +
    geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                      ymax = exp(y + 1.96*sqrt(pvar))),
                  width=.2) +
    labs(y="Species richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.fact.Species)

cowplot::plot_grid(plotlist=c(plots.cont.Species, plots.fact.Species))
ggsave("prediction.plots.Species.png", 
       plot = cowplot::plot_grid(plotlist=c(plots.cont.Species, plots.fact.Species)), 
       width = 27, height = 18, units = "cm", dpi=300)


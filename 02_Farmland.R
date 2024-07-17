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


# Farmland richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Farmland <- glmer(Farmland~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                     (1|Highway/transect), data=df, family="poisson")
vif(m)

# random effects selection
m2.Farmland <- update(m.Farmland, ~.-(1|Highway/transect)+(1|Highway))
m3.Farmland <- glm(Farmland~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc, data=df, family="poisson")
anova(m.Farmland, m2.Farmland, m3.Farmland)
m.Farmland <- glmer(Farmland~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                     (1|Highway), data=df, family="poisson")

# fixed effects selection
summary(m.Farmland)
m.Farmland <- stepAIC(m.Farmland)

# diagnostics
check_overdispersion(m.Farmland)

# model outputs -------------------
(s.Farmland <- summary(m.Farmland))
(drop.Farmland <- drop1(m.Farmland, test="Chisq"))
(r2.Farmland <- r2(m.Farmland))

(preds.Farmland <- all.vars(formula(m.Farmland))[-1][-6])
preds.cont.Farmland <- preds.Farmland[sapply(preds.Farmland, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Farmland <- preds.Farmland[sapply(preds.Farmland, function(pred) is.factor(pull(df, pred)))]
varimp.r2.Farmland <- sapply(preds.Farmland, function(pred){
  r2.Farmland$R2_marginal - r2(update(m.Farmland, paste("~.-", pred, sep="")))$R2_marginal
}) %>% set_names(preds.Farmland)
barplot(sort(varimp.r2.Farmland, decreasing = T), las=3)

# predictions
nds.cont.Farmland <- lapply(preds.cont.Farmland, function(pred){
  nd <- lapply(preds.cont.Farmland, function(pr) rep(mean(pull(df, pr)), 100)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Farmland) %>% 
    cbind(lapply(preds.fact.Farmland, function(pr) factor(rep(levels(pull(df, pr))[1], 100), levels=levels(pull(df, pr)))) %>% 
            bind_cols %>% 
            set_names(preds.fact.Farmland))
  nd[,pred] <- seq(min(pull(df, pred)),  max(pull(df, pred)), l=100)
  mm <- model.matrix(as.formula(paste("~", paste(preds.Farmland, collapse = "+"), sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Farmland)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Farmland), mm))
  nd
}) %>% set_names(preds.cont.Farmland)

nds.fact.Farmland <- lapply(preds.fact.Farmland, function(pred){
  nd <- lapply(preds.cont.Farmland, function(pr) rep(mean(pull(df, pr)), length(levels(pull(df, pred))))) %>% 
    bind_cols %>% 
    set_names(preds.cont.Farmland) %>% 
    cbind(lapply(preds.fact.Farmland, function(pr) 
      factor(rep(levels(pull(df, pr))[1], length(levels(pull(df, pred)))), levels=levels(pull(df, pr)))) %>% 
        bind_cols %>% 
        set_names(preds.fact.Farmland))
  nd[,pred] <- factor(levels(pull(df, pred)), levels = levels(pull(df, pred)))
  mm <- model.matrix(as.formula(paste("~", paste(preds.Farmland, collapse = "+"), sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Farmland)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Farmland), mm))
  nd
}) %>% set_names(preds.fact.Farmland)

# prediction plots
plots.cont.Farmland <- lapply(preds.cont.Farmland, function(pred){
  nd <- nds.cont.Farmland[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Farmland), alpha = .1) +
    geom_line() +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar))), alpha= .3) +
    labs(y="Farmland richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Farmland)
plots.fact.Farmland <- lapply(preds.fact.Farmland, function(pred){
  nd <- nds.fact.Farmland[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred))) +
    geom_jitter(data=df, aes(y=Farmland), alpha = .1, width=.2) +
    geom_point() +
    geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                      ymax = exp(y + 1.96*sqrt(pvar))),
                  width=.2) +
    labs(y="Farmland richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.fact.Farmland)

cowplot::plot_grid(plotlist=c(plots.cont.Farmland, plots.fact.Farmland))
ggsave("prediction.plots.Farmland.png", 
       plot = cowplot::plot_grid(plotlist=c(plots.cont.Farmland, plots.fact.Farmland)), 
       width = 27, height = 18, units = "cm", dpi=300)

outs.Farmland = list(
  model = m.Farmland,
  summary = s.Farmland,
  drop = drop.Farmland,
  r2 = r2.Farmland,
  varimp = varimp.r2.Farmland,
  preds = preds.Farmland,
  preds.cont = preds.cont.Farmland,
  preds.fact = preds.fact.Farmland,
  nds.cont = nds.cont.Farmland,
  nds.fact = nds.fact.Farmland
)
save(outs.Farmland, file = "outs.Farmland.RData")

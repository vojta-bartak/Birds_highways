library(tidyverse)
library(readxl)
library(sf)
library(lme4)
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


# Canopy richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Canopy <- glmer(Canopy~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                     (1|Highway/transect), data=df, family="poisson")
vif(m.Canopy)

# random effects selection
m2.Canopy <- update(m.Canopy, ~.-(1|Highway/transect)+(1|Highway))
m3.Canopy <- glm(Canopy~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc, data=df, family="poisson")
anova(m.Canopy, m2.Canopy, m3.Canopy)
m.Canopy <- glmer(Canopy~distance_sc+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                     (1|Highway), data=df, family="poisson")

# fixed effects selection
summary(m.Canopy)
m.Canopy <- stepAIC(m.Canopy)

# diagnostics
check_overdispersion(m.Canopy)

# model outputs -------------------
(s.Canopy <- summary(m.Canopy))
(drop.Canopy <- drop1(m.Canopy, test="Chisq"))
(r2.Canopy <- r2(m.Canopy))

(preds.Canopy <- all.vars(formula(m.Canopy))[-1][-6])
preds.cont.Canopy <- preds.Canopy[sapply(preds.Canopy, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Canopy <- preds.Canopy[sapply(preds.Canopy, function(pred) is.factor(pull(df, pred)))]
varimp.r2.Canopy <- sapply(preds.Canopy, function(pred){
  r2.Canopy$R2_marginal - r2(update(m.Canopy, paste("~.-", pred, sep="")))$R2_marginal
}) %>% set_names(preds.Canopy)
barplot(sort(varimp.r2.Canopy, decreasing = T), las=3)

# predictions
nds.cont.Canopy <- lapply(preds.cont.Canopy, function(pred){
  nd <- lapply(preds.cont.Canopy, function(pr) rep(mean(pull(df, pr)), 100)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Canopy) 
  nd[,pred] <- seq(min(pull(df, pred)),  max(pull(df, pred)), l=100)
  mm <- model.matrix(as.formula(paste("~", paste(preds.Canopy, collapse = "+"), sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Canopy)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Canopy), mm))
  nd
}) %>% set_names(preds.cont.Canopy)


# prediction plots
plots.cont.Canopy <- lapply(preds.cont.Canopy, function(pred){
  nd <- nds.cont.Canopy[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Canopy), alpha = .1) +
    geom_line() +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar))), alpha= .3) +
    labs(y="Canopy richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Canopy)

cowplot::plot_grid(plotlist=plots.cont.Canopy)
ggsave("prediction.plots.Canopy.png", 
       plot = cowplot::plot_grid(plotlist=plots.cont.Canopy), 
       width = 27, height = 18, units = "cm", dpi=300)

outs.Canopy = list(
  model = m.Canopy,
  summary = s.Canopy,
  drop = drop.Canopy,
  r2 = r2.Canopy,
  varimp = varimp.r2.Canopy,
  preds = preds.Canopy,
  preds.cont = preds.cont.Canopy,
  preds.fact = NULL,
  nds.cont = nds.cont.Canopy,
  nds.fact = NULL
)
save(outs.Canopy, file = "outs.Canopy.RData")

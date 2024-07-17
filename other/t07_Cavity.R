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
  mutate(forest_sc = scale(Forest) %>% as.numeric,
         village_sc = scale(Village) %>% as.numeric,
         distance_sc = scale(Distance) %>% as.numeric,
         traffic_sc = scale(Traffic_volume) %>% as.numeric,
         alt_sc = scale(Altitude) %>% as.numeric,
         e2_sc = scale(E2) %>% as.numeric,
         e3_sc = scale(E3) %>% as.numeric,
         water = as.factor(Water),
         ruderal = as.factor(Ruderal),
         meadow = as.factor(Meadow),
         field = as.factor(Field),
         noise = as.factor(Noise),
         Type = as.factor(Type),
         transect = sapply(ID, function(id) substr(id, 1, str_length(id)-1))
  )
orig.preds <- c(forest_sc="Forest", village_sc="Village", distance_sc="Distance", traffic_sc="Traffic_volume", alt_sc="Altitude",
                e2_sc="E2", e3_sc="E3", water="Water", ruderal="Ruderal", meadow="Meadow", field="Field", noise="Noise", Type="Type")

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

# Cavity richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Cavity.t <- glmer(Cavity~Type+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                        (1|Highway), data=df, family="poisson")
vif(m.Cavity.t)

m.Cavity.t <- glmer(Cavity~Type*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc)+
                        (1|Highway), data=df, family="poisson")


# random effects selection
m2.Cavity.t <- glm(Cavity~Type*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc), data=df, family="poisson")
anova(m.Cavity.t, m2.Cavity.t)

# fixed effects selection
summary(m.Cavity.t)
m.Cavity.t <- stepAIC(m2.Cavity.t)

# diagnostics
check_overdispersion(m.Cavity.t)

# model outputs -------------------
(s.Cavity <- summary(m.Cavity.t))
(drop.Cavity <- drop1(m.Cavity.t, test="Chisq"))
(r2.Cavity <- r2(m.Cavity.t))

preds.Cavity <- all.vars(formula(m.Cavity.t))[-1]
(preds.Cavity <- preds.Cavity[preds.Cavity!="Highway"])
preds.cont.Cavity <- preds.Cavity[sapply(preds.Cavity, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Cavity <- preds.Cavity[sapply(preds.Cavity, function(pred) is.factor(pull(df, pred)))]
varimp.r2.Cavity <- sapply(preds.Cavity, function(pred){
  r2.Cavity$R2_Nagelkerke - r2(update(m.Cavity.t, paste("~.-", pred,"-Type:",pred, sep="")))$R2_Nagelkerke
}) %>% set_names(preds.Cavity)
barplot(sort(varimp.r2.Cavity, decreasing = T), las=3)

# predictions
nds.cont.Cavity <- lapply(preds.cont.Cavity, function(pred){
  nd <- lapply(preds.cont.Cavity, function(pr) rep(mean(pull(df, pr)), 400)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Cavity) %>% 
    cbind(lapply(preds.fact.Cavity, function(pr) {
      factor(rep(levels(pull(df, pr))[1], 400), levels=levels(pull(df, pr)))
    }) %>% 
      bind_cols %>% 
      set_names(preds.fact.Cavity))
  nd[,pred] <- rep(seq(min(pull(df, pred)),  max(pull(df, pred)), l=100), 4)
  nd$Type <- factor(rep(unique(df$Type), each=100), levels=levels(df$Type))
  prediction <- predict(m.Cavity.t, newdata=nd, type="link", se=T)
  nd$y <- prediction$fit
  nd$pvar <- prediction$se.fit**2
  nd
}) %>% set_names(preds.cont.Cavity)

nds.fact.Cavity <- lapply(preds.fact.Cavity, function(pred){
  nd <- lapply(preds.cont.Cavity, function(pr) rep(mean(pull(df, pr)), length(levels(pull(df, pred))))) %>% 
    bind_cols %>% 
    set_names(preds.cont.Cavity) %>% 
    cbind(lapply(preds.fact.Cavity, function(pr) 
      factor(rep(levels(pull(df, pr))[1], length(levels(pull(df, pred)))), levels=levels(pull(df, pr)))) %>% 
        bind_cols %>% 
        set_names(preds.fact.Cavity))
  nd[,pred] <- factor(levels(pull(df, pred)), levels = levels(pull(df, pred)))
  prediction <- predict(m.Cavity.t, newdata=nd, type="link", se=T)
  nd$y <- prediction$fit
  nd$pvar <- prediction$se.fit**2
  nd
}) %>% set_names(preds.fact.Cavity)

# prediction plots
plots.cont.Cavity <- lapply(preds.cont.Cavity, function(pred){
  nd <- nds.cont.Cavity[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Cavity, color=Type), alpha = .1) +
    geom_line(aes(color=Type)) +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar)),
                    fill=Type), alpha= .3) +
    labs(y="Cavity richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Cavity)
plots.fact.Cavity <- lapply(preds.fact.Cavity, function(pred){
  nd <- nds.fact.Cavity[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred), color=Type)) +
    geom_jitter(data=df, aes(y=Cavity), alpha = .1, width=.2) +
    geom_point() +
    geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                      ymax = exp(y + 1.96*sqrt(pvar))),
                  width=.2) +
    labs(y="Cavity richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.fact.Cavity)

cowplot::plot_grid(plotlist=c(plots.cont.Cavity, plots.fact.Cavity))
ggsave("prediction.plots.Cavity.t.png", 
       plot = cowplot::plot_grid(plotlist=c(plots.cont.Cavity, plots.fact.Cavity)), 
       width = 27, height = 18, units = "cm", dpi=300)

outs.Cavity.t = list(
  model = m.Cavity.t,
  summary = s.Cavity,
  drop = drop.Cavity,
  r2 = r2.Cavity,
  varimp = varimp.r2.Cavity,
  preds = preds.Cavity,
  preds.cont = preds.cont.Cavity,
  preds.fact = preds.fact.Cavity,
  nds.cont = nds.cont.Cavity,
  nds.fact = nds.fact.Cavity
)
save(outs.Cavity.t, file = "t.outs.Cavity.RData")

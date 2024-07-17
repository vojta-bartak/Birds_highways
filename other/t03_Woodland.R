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

# Woodland richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Woodland.t <- glmer(Woodland~Type+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                        (1|Highway), data=df, family="poisson")
vif(m.Woodland.t)

m.Woodland.t <- glmer(Woodland~Type*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc)+
                        (1|Highway), data=df, family="poisson")


# random effects selection
m2.Woodland.t <- glm(Woodland~Type*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc), data=df, family="poisson")
anova(m.Woodland.t, m2.Woodland.t)

# fixed effects selection
summary(m.Woodland.t)
m.Woodland.t <- stepAIC(m2.Woodland.t)

# diagnostics
check_overdispersion(m.Woodland.t)

# model outputs -------------------
(s.Woodland <- summary(m.Woodland.t))
(drop.Woodland <- drop1(m.Woodland.t, test="Chisq"))
(r2.Woodland <- r2(m.Woodland.t))

preds.Woodland <- all.vars(formula(m.Woodland.t))[-1]
(preds.Woodland <- preds.Woodland[preds.Woodland!="Highway"])
preds.cont.Woodland <- preds.Woodland[sapply(preds.Woodland, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Woodland <- preds.Woodland[sapply(preds.Woodland, function(pred) is.factor(pull(df, pred)))]
varimp.r2.Woodland <- sapply(preds.Woodland, function(pred){
  r2.Woodland$R2_Nagelkerke - r2(update(m.Woodland.t, paste("~.-", pred,"-Type:",pred, sep="")))$R2_Nagelkerke
}) %>% set_names(preds.Woodland)
barplot(sort(varimp.r2.Woodland, decreasing = T), las=3)

# predictions
nds.cont.Woodland <- lapply(preds.cont.Woodland, function(pred){
  nd <- lapply(preds.cont.Woodland, function(pr) rep(mean(pull(df, pr)), 400)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Woodland) %>% 
    cbind(lapply(preds.fact.Woodland, function(pr) {
      factor(rep(levels(pull(df, pr))[1], 400), levels=levels(pull(df, pr)))
      }) %>% 
            bind_cols %>% 
            set_names(preds.fact.Woodland))
  nd[,pred] <- rep(seq(min(pull(df, pred)),  max(pull(df, pred)), l=100), 4)
  nd$Type <- factor(rep(unique(df$Type), each=100), levels=levels(df$Type))
  prediction <- predict(m.Woodland.t, newdata=nd, type="link", se=T)
  nd$y <- prediction$fit
  nd$pvar <- prediction$se.fit**2
  nd
}) %>% set_names(preds.cont.Woodland)

nds.fact.Woodland <- lapply(preds.fact.Woodland, function(pred){
  nd <- lapply(preds.cont.Woodland, function(pr) rep(mean(pull(df, pr)), length(levels(pull(df, pred))))) %>% 
    bind_cols %>% 
    set_names(preds.cont.Woodland) %>% 
    cbind(lapply(preds.fact.Woodland, function(pr) 
      factor(rep(levels(pull(df, pr))[1], length(levels(pull(df, pred)))), levels=levels(pull(df, pr)))) %>% 
        bind_cols %>% 
        set_names(preds.fact.Woodland))
  nd[,pred] <- factor(levels(pull(df, pred)), levels = levels(pull(df, pred)))
  prediction <- predict(m.Woodland.t, newdata=nd, type="link", se=T)
  nd$y <- prediction$fit
  nd$pvar <- prediction$se.fit**2
  nd
}) %>% set_names(preds.fact.Woodland)

# prediction plots
plots.cont.Woodland <- lapply(preds.cont.Woodland, function(pred){
  nd <- nds.cont.Woodland[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Woodland, color=Type), alpha = .1) +
    geom_line(aes(color=Type)) +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar)),
                    fill=Type), alpha= .3) +
    labs(y="Woodland richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Woodland)
plots.fact.Woodland <- lapply(preds.fact.Woodland, function(pred){
  nd <- nds.fact.Woodland[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred), color=Type)) +
    geom_jitter(data=df, aes(y=Woodland), alpha = .1, width=.2) +
    geom_point() +
    geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                      ymax = exp(y + 1.96*sqrt(pvar))),
                  width=.2) +
    labs(y="Woodland richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.fact.Woodland)

cowplot::plot_grid(plotlist=c(plots.cont.Woodland, plots.fact.Woodland))
ggsave("prediction.plots.Woodland.t.png", 
       plot = cowplot::plot_grid(plotlist=c(plots.cont.Woodland, plots.fact.Woodland)), 
       width = 27, height = 18, units = "cm", dpi=300)

outs.Woodland.t = list(
  model = m.Woodland.t,
  summary = s.Woodland,
  drop = drop.Woodland,
  r2 = r2.Woodland,
  varimp = varimp.r2.Woodland,
  preds = preds.Woodland,
  preds.cont = preds.cont.Woodland,
  preds.fact = preds.fact.Woodland,
  nds.cont = nds.cont.Woodland,
  nds.fact = nds.fact.Woodland
)
save(outs.Woodland.t, file = "t.outs.Woodland.RData")

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

# Species richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Species.t <- glmer(Species~Type+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                     (1|Highway), data=df, family="poisson")
vif(m.Species.t)

m.Species.t <- glmer(Species~Type*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc)+
                       (1|Highway), data=df, family="poisson")


# random effects selection
m2.Species.t <- glm(Species~Type*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc), data=df, family="poisson")
anova(m.Species.t, m2.Species.t)

# fixed effects selection
summary(m.Species.t)
m.Species.t <- stepAIC(m.Species.t)

# diagnostics
check_overdispersion(m.Species.t)

# model outputs -------------------
(s.Species <- summary(m.Species.t))
(drop.Species <- drop1(m.Species.t, test="Chisq"))
(r2.Species <- r2(m.Species.t))

(preds.Species <- all.vars(formula(m.Species.t))[-1][-10])
preds.cont.Species <- preds.Species[sapply(preds.Species, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Species <- preds.Species[sapply(preds.Species, function(pred) is.factor(pull(df, pred)))]
varimp.r2.Species <- sapply(preds.Species, function(pred){
  r2.Species$R2_marginal - r2(update(m.Species.t, paste("~.-", pred,"-Type:",pred, sep="")))$R2_marginal
}) %>% set_names(preds.Species)
barplot(sort(varimp.r2.Species, decreasing = T), las=3)

# predictions
nds.cont.Species <- lapply(preds.cont.Species, function(pred){
  nd <- lapply(preds.cont.Species, function(pr) rep(mean(pull(df, pr)), 400)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Species) %>% 
    cbind(lapply(preds.fact.Species, function(pr) factor(rep(levels(pull(df, pr))[1], 400), levels=levels(pull(df, pr)))) %>% 
            bind_cols %>% 
            set_names(preds.fact.Species))
  nd[,pred] <- rep(seq(min(pull(df, pred)),  max(pull(df, pred)), l=100), 4)
  nd$Type <- factor(rep(unique(df$Type), each=100), levels=levels(df$Type))
  mm <- model.matrix(as.formula(paste("~", 
                                      paste(
                                        as.character(formula(m.Species.t))[3] %>% 
                                        strsplit(split = " + ", fixed=T) %>% 
                                        "[["(1) %>% 
                                        "["(. != "(1 | Highway)"),
                                        collapse = "+"), 
                                      sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Species.t)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Species.t), mm))
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
  mm <- model.matrix(as.formula(paste("~", 
                                      paste(
                                        as.character(formula(m.Species.t))[3] %>% 
                                          strsplit(split = " + ", fixed=T) %>% 
                                          "[["(1) %>% 
                                          "["(. != "(1 | Highway)"),
                                        collapse = "+"), 
                                      sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Species.t)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Species.t), mm))
  nd
}) %>% set_names(preds.fact.Species)

# prediction plots
plots.cont.Species <- lapply(preds.cont.Species, function(pred){
  nd <- nds.cont.Species[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Species, color=Type), alpha = .1) +
    geom_line(aes(color=Type)) +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar)),
                    fill=Type), alpha= .3) +
    labs(y="Species richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Species)
plots.fact.Species <- lapply(preds.fact.Species, function(pred){
  nd <- nds.fact.Species[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred), color=Type)) +
    geom_jitter(data=df, aes(y=Species), alpha = .1, width=.2) +
    geom_point() +
    geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                      ymax = exp(y + 1.96*sqrt(pvar))),
                  width=.2) +
    labs(y="Species richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.fact.Species)

cowplot::plot_grid(plotlist=c(plots.cont.Species, plots.fact.Species))
ggsave("prediction.plots.Species.t.png", 
       plot = cowplot::plot_grid(plotlist=c(plots.cont.Species, plots.fact.Species)), 
       width = 27, height = 18, units = "cm", dpi=300)

outs.Species.t = list(
  model = m.Species.t,
  summary = s.Species,
  drop = drop.Species,
  r2 = r2.Species,
  varimp = varimp.r2.Species,
  preds = preds.Species,
  preds.cont = preds.cont.Species,
  preds.fact = preds.fact.Species,
  nds.cont = nds.cont.Species,
  nds.fact = nds.fact.Species
)
save(outs.Species.t, file = "t.outs.Species.RData")

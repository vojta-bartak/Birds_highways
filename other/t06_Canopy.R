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

# Canopy richness -------------------------------------------------------------------------------------------------------------

# full model (with traffic volume)
m.Canopy.t <- glmer(Canopy~Type+e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc+
                        (1|Highway), data=df, family="poisson")
vif(m.Canopy.t)

m.Canopy.t <- glmer(Canopy~Type*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc)+
                        (1|Highway), data=df, family="poisson")


# random effects selection
m2.Canopy.t <- glm(Canopy~Type*(e2_sc+e3_sc+field+meadow+ruderal+water+forest_sc+village_sc+traffic_sc+alt_sc), data=df, family="poisson")
anova(m.Canopy.t, m2.Canopy.t)

# fixed effects selection
summary(m.Canopy.t)
m.Canopy.t <- stepAIC(m.Canopy.t)

# diagnostics
check_overdispersion(m.Canopy.t)

# model outputs -------------------
(s.Canopy <- summary(m.Canopy.t))
(drop.Canopy <- drop1(m.Canopy.t, test="Chisq"))
(r2.Canopy <- r2(m.Canopy.t))

preds.Canopy <- all.vars(formula(m.Canopy.t))[-1]
(preds.Canopy <- preds.Canopy[preds.Canopy!="Highway"])
preds.cont.Canopy <- preds.Canopy[sapply(preds.Canopy, function(pred) is.numeric(pull(df, pred)))]
preds.fact.Canopy <- preds.Canopy[sapply(preds.Canopy, function(pred) is.factor(pull(df, pred)))]
varimp.r2.Canopy <- sapply(preds.Canopy, function(pred){
  r2.Canopy$R2_marginal - r2(update(m.Canopy.t, paste("~.-", pred,"-Type:",pred, sep="")))$R2_marginal
}) %>% set_names(preds.Canopy)
barplot(sort(varimp.r2.Canopy, decreasing = T), las=3)

# predictions
nds.cont.Canopy <- lapply(preds.cont.Canopy, function(pred){
  nd <- lapply(preds.cont.Canopy, function(pr) rep(mean(pull(df, pr)), 400)) %>% 
    bind_cols %>% 
    set_names(preds.cont.Canopy) %>% 
    cbind(lapply(preds.fact.Canopy, function(pr) factor(rep(levels(pull(df, pr))[1], 400), levels=levels(pull(df, pr)))) %>% 
            bind_cols %>% 
            set_names(preds.fact.Canopy))
  nd[,pred] <- rep(seq(min(pull(df, pred)),  max(pull(df, pred)), l=100), 4)
  nd$Type <- factor(rep(unique(df$Type), each=100), levels=levels(df$Type))
  mm <- model.matrix(as.formula(paste("~", 
                                      paste(
                                        as.character(formula(m.Canopy.t))[3] %>% 
                                          strsplit(split = " + ", fixed=T) %>% 
                                          "[["(1) %>% 
                                          "["(. != "(1 | Highway)"),
                                        collapse = "+"), 
                                      sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Canopy.t)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Canopy.t), mm))
  nd
}) %>% set_names(preds.cont.Canopy)

nds.fact.Canopy <- lapply(preds.fact.Canopy, function(pred){
  nd <- lapply(preds.cont.Canopy, function(pr) rep(mean(pull(df, pr)), length(levels(pull(df, pred))))) %>% 
    bind_cols %>% 
    set_names(preds.cont.Canopy) %>% 
    cbind(lapply(preds.fact.Canopy, function(pr) 
      factor(rep(levels(pull(df, pr))[1], length(levels(pull(df, pred)))), levels=levels(pull(df, pr)))) %>% 
        bind_cols %>% 
        set_names(preds.fact.Canopy))
  nd[,pred] <- factor(levels(pull(df, pred)), levels = levels(pull(df, pred)))
  mm <- model.matrix(as.formula(paste("~", 
                                      paste(
                                        as.character(formula(m.Canopy.t))[3] %>% 
                                          strsplit(split = " + ", fixed=T) %>% 
                                          "[["(1) %>% 
                                          "["(. != "(1 | Highway)"),
                                        collapse = "+"), 
                                      sep="")), data=nd)
  nd$y <- mm %*% fixef(m.Canopy.t)
  nd$pvar <- diag(mm %*% tcrossprod(vcov(m.Canopy.t), mm))
  nd
}) %>% set_names(preds.fact.Canopy)

# prediction plots
plots.cont.Canopy <- lapply(preds.cont.Canopy, function(pred){
  nd <- nds.cont.Canopy[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred)*sd(pull(df, orig.preds[pred])) + mean(pull(df, orig.preds[pred])))) +
    geom_point(data=df, aes(y=Canopy, color=Type), alpha = .1) +
    geom_line(aes(color=Type)) +
    geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                    ymax = exp(y + 1.96*sqrt(pvar)),
                    fill=Type), alpha= .3) +
    labs(y="Canopy richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.cont.Canopy)
plots.fact.Canopy <- lapply(preds.fact.Canopy, function(pred){
  nd <- nds.fact.Canopy[[pred]]
  ggplot(nd, aes(y=exp(y), x=get(pred), color=Type)) +
    geom_jitter(data=df, aes(y=Canopy), alpha = .1, width=.2) +
    geom_point() +
    geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)),
                      ymax = exp(y + 1.96*sqrt(pvar))),
                  width=.2) +
    labs(y="Canopy richness", x=orig.preds[pred]) +
    theme_bw()
}) %>% set_names(preds.fact.Canopy)

cowplot::plot_grid(plotlist=c(plots.cont.Canopy, plots.fact.Canopy))
ggsave("prediction.plots.Canopy.t.png", 
       plot = cowplot::plot_grid(plotlist=c(plots.cont.Canopy, plots.fact.Canopy)), 
       width = 27, height = 18, units = "cm", dpi=300)

outs.Canopy.t = list(
  model = m.Canopy.t,
  summary = s.Canopy,
  drop = drop.Canopy,
  r2 = r2.Canopy,
  varimp = varimp.r2.Canopy,
  preds = preds.Canopy,
  preds.cont = preds.cont.Canopy,
  preds.fact = preds.fact.Canopy,
  nds.cont = nds.cont.Canopy,
  nds.fact = nds.fact.Canopy
)
save(outs.Canopy.t, file = "t.outs.Canopy.RData")

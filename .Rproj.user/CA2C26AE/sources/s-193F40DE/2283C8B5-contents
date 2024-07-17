library(tidyverse)
library(readxl)

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
orig.preds <- c(distance_sc="Distance", e3_sc="E3", e2_sc="E2", traffic_sc="Traffic_volume", alt_sc="Altitude",
                forest_sc="Forest", village_sc="Village", water="Water", ruderal="Ruderal", meadow="Meadow", field="Field")
label_preds <- c("Dist. to motorway (m)", "Tree cover (%)", "Shrub cover (%)", "Traffic volume (cars per day)", "Elevation (m)",
                 "Dist. to forest (m)", "Dist. to village (m)", "Water", "Ruderal", "Meadow", "Arrable")

groups <- data.frame(response = c("Species","Farmland", "Woodland", "Ground", "Shrub", "Canopy", 
                                  "Cavity", "Insect", "Granivor"),
                     type = c("All", "Landsclape", "Landscape", "Nesting", "Nesting", "Nesting", 
                              "Nesting", "Diet", "Diet"),
                     label = c("All species", "Farmland", "Woodland", "Ground", "Shrub", "Canopy", 
                               "Cavity", "Insectivores", "Granivores"))
for (resp in groups$response) load(paste("outs",resp,"RData", sep="."))

# table of coefficients -------------------------------------------------------------------------------

# anova tables (significances) ------------------------------------------------------------------------
drops <- lapply(groups$response, function(resp){
  get(paste("outs",resp,sep="."))$drop %>% 
    mutate(response = resp,
           p = if ("Pr(>Chi)" %in% names(.)) get("Pr(>Chi)") else get("Pr(Chi)"),
           sign = p < 0.05) %>% 
    set_names(paste(names(.), pull(., "response")[1], sep=".")) %>% 
    rownames_to_column("predictor")
}) %>% 
  reduce(full_join, by="predictor") %>% 
  filter(predictor != "<none>")

drops %>% select(predictor, contains("sign")) %>% write.table("anovas.1.csv", row.names = F, sep=";")


# variable importances --------------------------------------------------------------------------------
importances <- lapply(1:nrow(groups), function(i){
  resp <- groups$response[i]
  get(paste("outs",resp,sep="."))$varimp %>% 
    as.data.frame() %>% 
    set_names(resp) %>% 
    rownames_to_column("predictor")
}) %>% 
  reduce(left_join, by="predictor") %>% 
  pivot_longer(2:10) %>%
  mutate(name = factor(name, levels = groups$response, labels = groups$label),
         value = ifelse(value > 0, value, 0),
         predictor = factor(predictor, 
                            levels = c("distance_sc", "traffic_sc", "alt_sc","e3_sc", "e2_sc", 
                                       "village_sc", "field", "water", "ruderal"),
                            labels = c("Dist. to motorway", "Traffic volume", "Altitude", "Tree cover", "Shrub cover", 
                                       "Dist. to village", "Arrable", "Water", "Ruderal")))
  
ggplot(importances, aes(x=predictor, y=value*100)) +
  geom_col() +
  facet_wrap(~name) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=.5)) +
  labs(y = "Explained variability (%)", x="")
ggsave("importances.png", dpi=300, width=18, height=16, units = "cm")



# r2 values -------------------------------------------------------------------------------------------
r2s <- sapply(groups$response, function(resp){
  rr = get(paste("outs",resp,sep="."))$r2
  rr[length(rr)][[1]]
}) %>% set_names(groups$response) 

r2s %>% barplot(las = 3)
types <- sapply(groups$response, function(resp){
  rr = get(paste("outs",resp,sep="."))$r2
  if (names(rr)[1] == "R2_conditional") "GLMM" else "GLM"
}) %>% set_names(groups$response) 

models <- lapply(groups$response, function(resp){
  get(paste("outs",resp,sep="."))$model
})
m <- models[[1]]
slots(m)
strc <- sapply(models, function(m) formula(m)[[3]] %>% as.character %>% strsplit(split=" + ", fixed=T) %>% "[["(3))
sfits <- c(Species=F, Farmland=F, Woodland=F, Ground=T, Shrub=F, Canopy=T, Cavity=T, Insect=F, Granivor=T)

# model summaries -------------------------------------------------------------------------------------
data.frame(model = names(r2s), type = types, r2 = round(r2s*100, 1), ranef = ifelse(strc!="alt_sc",strc,"-"), singular = sfits) %>% 
  write.table("model_summaries.csv", sep=";", row.names = F)

# table of coefficients -------------------------------------------------------------------------------
coefs <- lapply(groups$response, function(resp){
  s <- get(paste("outs",resp,sep="."))$summary$coefficients %>% as.data.frame
  alpha <- case_when(s[,"Pr(>|z|)"] < 0.001 ~ "***",
                     s[,"Pr(>|z|)"] < 0.01 ~ "**",
                     s[,"Pr(>|z|)"] < 0.05 ~ "*",
                     s[,"Pr(>|z|)"] < 0.1 ~ "\u00B7",
                     TRUE ~ "")
  s[,resp] <- paste(round(s[,"Estimate"], 3), "\u00B1", round(s[,"Std. Error"], 3), alpha, sep="")
  s %>% rownames_to_column("predictor") %>% "["(,c("predictor", resp))
}) %>% reduce(full_join)
coefs[is.na(coefs)] <- "-"
write.table(coefs, "coefficients.csv", sep=";", row.names = F)

# prediction plots ------------------------------------------------------------------------------------
nds.cont <- lapply(groups$response, function(resp){
  nds <- get(paste("outs",resp,sep="."))$nds.cont
  lapply(names(nds), function(pred){
    if (drops[drops$predictor==pred, paste("sign",resp,sep=".")]){
      nds[[pred]] %>% 
        mutate(x = get(pred)*sd(pull(df, orig.preds[pred]))+mean(pull(df, orig.preds[pred])),
               pred = pred,
               resp = resp)
    } else {
      NULL
    }
      
  }) %>% 
    bind_rows
}) %>% bind_rows

nds.fact <- lapply(groups$response, function(resp){
  nds <- get(paste("outs",resp,sep="."))$nds.fact
  lapply(names(nds), function(pred){
    if (drops[drops$predictor==pred, paste("sign",resp,sep=".")]){
      nds[[pred]] %>% 
        mutate(x = get(pred),
               pred = pred,
               resp = resp)
    } else {
      NULL
    }
    
  }) %>% 
    bind_rows
}) %>% bind_rows

# All groups
nds.cont %>% 
  ggplot(aes(x = x, y = exp(y), color = resp, fill = resp)) +
  geom_line() +
  geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)), ymax = exp(y + 1.96*sqrt(pvar))), alpha=.3) +
  facet_wrap(~pred, scales = "free_x")
nds.fact %>% 
  ggplot(aes(x = x, y = exp(y), color = resp)) +
  geom_point(position = position_dodge(width=.2)) +
  geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)), ymax = exp(y + 1.96*sqrt(pvar))), position = position_dodge(width = .2)) +
  facet_wrap(~pred)

# Farmland / Woodland
p1.Farm.Wood <- nds.cont %>% 
  filter(resp %in% c("Species", "Farmland", "Woodland")) %>% 
  mutate(resp = factor(resp, levels = c("Species", "Farmland", "Woodland"), labels = c("All species", "Farmland", "Woodland")),
         pred = factor(pred, levels = names(orig.preds), labels = label_preds)) %>% 
  ggplot(aes(x = x, y = exp(y))) +
  geom_line(aes(color = resp)) +
  geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)), ymax = exp(y + 1.96*sqrt(pvar)), fill = resp), alpha=.3) +
  facet_wrap(~pred, scales = "free_x") +
  labs(y = "Species richness", color="", fill="", x="")
ggsave("predictions.cont.Farm.Wood.png", height = 12, width = 16, units = "cm", dpi=300)
p2.Farm.Wood <- nds.fact %>% 
  filter(resp %in% c("Species", "Farmland", "Woodland")) %>% 
  mutate(resp = factor(resp, levels = c("Species", "Farmland", "Woodland"), labels = c("All species", "Farmland", "Woodland")),
         pred = factor(pred, levels = names(orig.preds), labels = orig.preds)) %>% 
  ggplot(aes(x = factor(x, levels=c(1,0), labels=c("Present", "Not present")), y = exp(y), color = resp)) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)), ymax = exp(y + 1.96*sqrt(pvar))), width=.2) +
  facet_wrap(~pred) +
  labs(y = "Species richness", color="", fill="", x="")
ggsave("predictions.fact.Farm.Wood.png", height = 12, width = 16, units = "cm", dpi=300)

cowplot::plot_grid(p1.Farm.Wood, cowplot::plot_grid(p2.Farm.Wood, NULL, nrow=1, rel_widths = c(2.8,1)), nrow=2, rel_heights = c(2,1))
ggsave("prediction_plots.Farm.Wood.png", height = 18, width = 18, units = "cm", dpi=300)

# Canopy / Cavity / Ground / Shrub
group <- c("Species", "Canopy", "Cavity", "Ground", "Shrub")
p1.Nest <- nds.cont %>% 
  filter(resp %in% group) %>% 
  mutate(resp = factor(resp, levels = group, labels = c("All species", group[2:length(group)])),
         pred = factor(pred, levels = names(orig.preds), labels = label_preds)) %>% 
  ggplot(aes(x = x, y = exp(y))) +
  geom_line(aes(color = resp)) +
  geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)), ymax = exp(y + 1.96*sqrt(pvar)), fill = resp), alpha=.3) +
  facet_wrap(~pred, scales = "free_x") +
  labs(y = "Species richness", color="", fill="", x="")
ggsave("predictions.cont.Nesting.png", height = 12, width = 16, units = "cm", dpi=300)
p2.Nest <- nds.fact %>% 
  filter(resp %in% group) %>% 
  mutate(resp = factor(resp, levels = group, labels = c("All species", group[2:length(group)])),
         pred = factor(pred, levels = names(orig.preds), labels = orig.preds)) %>% 
  ggplot(aes(x = factor(x, levels=c(1,0), labels=c("Present", "Not present")), y = exp(y), color = resp)) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)), ymax = exp(y + 1.96*sqrt(pvar))), width=.2) +
  facet_wrap(~pred) +
  labs(y = "Species richness", color="", fill="", x="")
ggsave("predictions.fact.Nesting.png", height = 12, width = 16, units = "cm", dpi=300)

cowplot::plot_grid(p1.Nest, p2.Nest, nrow=2, rel_heights = c(2,1))
ggsave("prediction_plots.Nest.png", height = 18, width = 18, units = "cm", dpi=300)

# Granivor / Insect
group <- c("Species", "Granivor", "Insect")
p1.Diet <- nds.cont %>% 
  filter(resp %in% group) %>% 
  mutate(resp = factor(resp, levels = group, labels = c("All species", "Granivors", "Insectivors")),
         pred = factor(pred, levels = names(orig.preds), labels = label_preds)) %>% 
  ggplot(aes(x = x, y = exp(y))) +
  geom_line(aes(color = resp)) +
  geom_ribbon(aes(ymin = exp(y - 1.96*sqrt(pvar)), ymax = exp(y + 1.96*sqrt(pvar)), fill = resp), alpha=.3) +
  facet_wrap(~pred, scales = "free_x") +
  labs(y = "Species richness", color="", fill="", x="")
ggsave("predictions.cont.Diet.png", height = 12, width = 16, units = "cm", dpi=300)
p2.Diet <- nds.fact %>% 
  filter(resp %in% group) %>% 
  mutate(resp = factor(resp, levels = group, labels = c("All species", "Granivors", "Insectivors")),
         pred = factor(pred, levels = names(orig.preds), labels = orig.preds)) %>% 
  ggplot(aes(x = factor(x, levels=c(1,0), labels=c("Present", "Not present")), y = exp(y), color = resp)) +
  geom_point() +
  geom_errorbar(aes(ymin = exp(y - 1.96*sqrt(pvar)), ymax = exp(y + 1.96*sqrt(pvar))), width=.2) +
  facet_wrap(~pred) +
  labs(y = "Species richness", color="", fill="", x="")
ggsave("predictions.fact.Diet.png", height = 12, width = 16, units = "cm", dpi=300)

cowplot::plot_grid(p1.Diet, cowplot::plot_grid(p2.Diet, NULL, nrow=1, rel_widths = c(2.8, 1)), nrow=2, rel_heights = c(2,1))
ggsave("prediction_plots.Diet.png", height = 18, width = 18, units = "cm", dpi=300)

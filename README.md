# Birds_highways

This repository contains a code and data for the following paper:

**Bird richness loss around motorways in fragmented agricultural landscape: a guild approach**

Hladík Štěpán a, ⁎, Zasadil Petr b, Barták Vojtěch c, Keken Zdeněk a

a *Department of Applied Ecology, Faculty of Environmental Sciences, Czech University of Life Sciences, Kamýcká 129, 165 00 Prague - Suchdol, Czech Republic*

*b* *Department of Ecology, Faculty of Environmental Sciences, Czech University of Life Sciences, Kamýcká 129, 165 00 Prague - Suchdol, Czech Republic*

*c* *Department of Spatial Sciences, Faculty of Environmental Sciences, Czech University of Life Sciences, Kamýcká 129, 165 00 Prague - Suchdol, Czech Republic*

⁎ Corresponding author.



The data are located in the file "data.xlsx". The code is organized as follows:

1. In each of the nine files "01_Species.R", "02_Farmland.R", ..., "09_Granivor.R", an analysis is made for a corresponding group of species. In each analysis, the data are loaded, models are fitted, model selection and diagnostics are performed, and model summaries are computed. The model outputs are then saved into a corresponding ".RData" file, i.e. "outs.Species.RData", "outs.Farmland.RData", etc.
2. In the file "outputs.R", the graphical and tabular outputs are made based on the model outputs saved in the previous step.
3. In the file "model_strategies.R", the competing model structures are compared based on AIC. The model strategies primarily differ in the way the traffic noise is represented, i.e. either by a distance to the nearest highway (predictor "distance_sc") or by discrete noise-level categories (predictor "noise"). See the main text of the paper for more details.
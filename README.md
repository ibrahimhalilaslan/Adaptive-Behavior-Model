This repository contains all the R code necessary to simulate the dynamics of schistosomiasis using a model that incorporates seasonality, with and without aestivation. The simulations are designed to evaluate the impact of seasonality and adaptive behavior on schistosomiasis and also predict the prevalence of S. haematobium and S. mansoni across Africa.

The repository is organized into two main folders:

1️⃣ aestivation/

Contains the simulation models with aestivation. It includes separate subfolders for S. haematobium and S. mansoni, each with the following structure:

* s.haematobium/
	
 •	🗂️ estimate_aestivation_parameters/
	•	R scripts for estimating aestivation function parameters.
	•	Includes maximum likelihood estimation (MLE) and AIC calculations for model selection.
	•	🗂️ heat_map/
	•	Scripts to generate heat maps for visualizing simulation results with GNTD data.
	•	Includes data sources used for the MLE function.
	•	🗂️ projection/
	•	Scripts to simulate and map schistosomiasis prevalence across Africa.
	•	Generates visual projections using GNTD data.
	•	🗂️ multiple_epsilon_simulation/
	•	R code for simulating the thermal performance curve with three different values of seasonality (epsilon).
	•	Helps analyze the impact of seasonal variation on schistosomiasis dynamics.
	
 
 * s.mansoni/

 •	Contains the same folder structure and functionality as s.haematobium/, but for S. mansoni simulations.

2️⃣ no_aestivation/

Contains the simulation models without aestivation. It also has separate subfolders for S. mansoni and S. haematobium, but with a simpler structure:

* s.haematobium/
	
 •	🗂️ heat_map/
	•	Scripts for generating heat maps of the model with seasonality, but without aestivation.
	•	Includes visual comparisons of simulation results with GNTD data.
	•	🗂️ model_run_with_different_seasonality/
	•	R scripts to simulate schistosomiasis prevalence using the model with seasonality (without aestivation).
	•	Simulations are run for three different seasonality values to observe the effects of varying seasonal intensity.
	
 *	s.mansoni/
   
	•	Contains the same folder structure and functionality as s.haematobium/, but for S. mansoni simulations.

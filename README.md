# HerbivoreDetritivoreInteractions
A repository for experimental data and code associated with the manuscript: Interactions between herbivores and detritivores develop slowly and conspicuously into feedbacks

The scripts in this repository replicate the statistical and theoretical analysis presented in the above manuscript. The order of analysis follows the order in the manuscript whenever possible.

1. The statistical models are all run in the script Statistical_Analysis.R. Section 1 is the analysis of the experiment. Section 2 is the analysis of the 1-m^2 plots that are reported in the appendicies.

2. The simple model and the comparison with the complex model using model results created elsewhere are run in simple_model.R

3. The complex model is run to equilibrium for versions with and without herbivores and detritivores in the script complex_model_eqm_to_cluster.R

4. The complex model is simulated based on empirical treatments so it can be compared to the data using the script complex_model_non-eqm_to_cluster.R

5. The analysis of the complex model from the script run in #4 and its fit to the empirical data is conducted in analysis_complex_model_post-cluster.R

6. The complex model with the subversions "direct model" and "multiple species model" are run in this code. Only one parameter set is run for each model in this code.

The other files contain functions that prepare raw data for analysis (extract_data_model_analysis.R), analyze climate data (extract_climate_data.R), and provide a convenient plotting function (plotcompare.R).

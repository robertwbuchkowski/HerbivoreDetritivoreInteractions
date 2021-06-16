# HerbivoreDetritivoreInteractions
A repository for experimental data and code associated with the manuscript: Interactions between herbivores and detritivores develop slowly and conspicuously into feedbacks

The scripts in this repository replicate the statistical and theoretical analysis presented in the above manuscript. The order of analysis follows the order in the manuscript whenever possible.

1. The statistical models are all run in the script Statistical_Analysis.R. Section 1 is the analysis of the experiment. Section 2 is the analysis of the 1-m^2 plots that are reported in the appendicies.

2. The model is simulated based on empirical treatments so it can be compared to the data using the script model_simulations_cluster.R

3. The analysis of the model from the script run in #2 and its fit to the empirical data is conducted in model_cluster_fit.R

4. The analysis of the model from the script run in #3 is in field_simulation_interaction_comparison.R. This analysis produces the data for the comparison with the field simulation and the plot of interaction strength in the best fitting model result.

5. The script model_simulations.R has the comparison between reverse and forward MM models and the long-term model simulations.

The other files contain functions that prepare raw data for analysis (extract_data_model_analysis.R), analyze climate data (extract_climate_data.R), and provide a convenient plotting function (plotcompare.R).

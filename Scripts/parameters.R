# Load in the parameters and staring values for the complex model:


# Baseline parameter values
params<- c(Vlm_mod = 8e-6,
           Vsm_mod = 4e-06,
           Klm_mod = 0.143,
           Ksm_mod = 0.143,
           Vlw = 2.4e-06,#2.4e-06, #2.4e-05 before correction to Type I
           Vsw = 4.1e-05,#4.1e-05,#0.00462 before correction to Type I
           Vpf = 0.001/0.00018, #From JRS project 0.03
           Kpf = 0.08, #From JRS project 0.006
           Vhp = 0.01, #From JRS project 0.0025 - 0.0029
           SUEh = 0.7,
           SUE = 0.50,
           SUEws = 0.01,
           SUEwl = 0.02,
           SUEwm = 0.3,
           q = 0.1,
           IN= 0.02,
           l = 0.0001,
           tm = 0.01,
           tw = 0.00001,
           th = 20, # based on surivival from Schmitz lab experiments
           fi=0.6, #0.6,
           fo=0.001, #0.003,
           tp = 0.000005, #0.00008,
           Ea = 0.25,
           Kappa = 8.62e-05,
           B = 2493,
           D = 26712,
           Vslope = 0.063,
           Kslope = 0.007,
           Vint = 5.47,
           Kint = 3.19,
           Tref_W = 288.15,
           Tref_P = 297.65)

# Starting inputs--> WE plots are the 1-m^2 plots discussed in the text
yint= c(P=95.238095, # WE with R:S ratio from Buchkowski et al. 2018
        L=25.526097, # WE plots with C:N ratio from files
        M=7.160995, # WE plots
        W=9.639477, # WE plots
        N=0.100000, # WE plots
        S=134.845515, # my historical data
        H=0.009) # Schmitz et al. 1997 8-10/m2 * 0.0986 * 0.11

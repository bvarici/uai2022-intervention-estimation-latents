
# intervention-estimation-latents
Codes for UAI 2022 paper: Intervention Target Estimation in the Presence of Latent Variables

**main_codes**: contain the code-base for the proposed PreDITEr algorithm.

run_simulations.py: The simulations with synthetic data. 
run_sachs_data.py: The simulations with protein signaling data (Sachs dataset). 
Preprocessed data is taken from https://github.com/csquires/utigsp.

plot_simulations.py: Generate Figure 2 in main text and Figure 4 in Appendix D.1.
plot_sachs_results.py: Generate Figure 3 in main text and Figure 7 in Appendix D.3.
(Requires graphviz package). 

Note: requires causaldag package (https://github.com/uhlerlab/causaldag).

TO-DO: create a notebook for step by step instructions.

Instructions for running the comparisons to $\psi$-FCI (Jaber et al. (2020)) and FCI-JCI123 (Mooij et al.(2020)) are given in "simulate_example_small_graph" (for a 4-node graph), and "simulate_larger_graphs" (for larger graphs) folders.

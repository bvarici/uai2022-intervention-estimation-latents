Compare our algorithm, psi-FCI and FCI-JCI123 methods on the example small graph.

*running our algorithm*

my_generate_data_cmd.py: generates gaussian data for the small graph.
my_execute_generate_data.sh: calls above function.
my_discretize_gaussian_data.py: discretize the gaussian data generated above.
my_run_discretized_gaussian.py: run our algorithm on the discretized gaussian data from the small graph.
my_plot_graphs: plot the graphs for results of our algorithm.


*running psi-FCI of Jaber et al.(2020) and FCI-JCI123 of Mooij et al. (2020)*

the auxillary codes for running these algorithms are obtained from the authors of psi-FCI paper.

run execute_psiFCI_gauss_discretized.sh and execute_JCI_gauss_discretized.sh files

All results are saved to folders with corresponding names.
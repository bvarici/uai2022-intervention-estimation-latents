Compare our algorithm and FCI-JCI123 methods on larger graphs.

the codes for running JCI algorithm of Mooij et al. (2020) is obtained from https://github.com/caus-am/jci.

History of running commands.

# generate data for JCI setting.
sh jci_generate_data.sh 10 1 0.2 0.2 5000 1 0 11 40
sh jci_generate_data.sh 20 1 0.1 0.1 5000 1 0 11 40
sh jci_generate_data.sh 30 1 0.066 0.066 5000 1 0 11 40
sh jci_generate_data.sh 40 1 0.05 0.05 5000 1 0 11 40 

execute "sh execute_my_generate_data.sh" with the following settings.

# jci-more-latents data generation.
python my_generate_data.py -r $run -s $sample -pSystem $pSystem -pContext $pContext -l 10 -i 1 -c 2

# jci-more-latents-more-targets data generation. (most logical)
python generate_data_flex.py -r $run -s $sample -pSystem $pSystem -pContext $pContext -l $pSystem -i 3 -c 2

# jci-multiple-contexts data generation.
python generate_data_flex.py -r $run -s $sample -pSystem $pSystem -pContext $pContext -l 5 -i 1 -c 2


*running our algorithm*
analyze_ours.py: imports necessary functions and runs our algorithm.

*running JCI*
under  "../experiments/simul/",

jci_execute_original.sh: runs for the JCI-generated data.
jci_execute_mydata-more-latents-more-targets.sh: runs for our generated data.


All results are saved to folders with corresponding names.
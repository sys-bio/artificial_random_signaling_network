# aritificial_random_signaling_network
If you use any of the codes, please cite the bioRxiv (https://www.biorxiv.org/content/10.1101/2020.05.08.084848v1) and the GitHub website (https://github.com/SunnyXu/artificial_random_signaling_network).

This package was implemented in Julia 1.2 on Windows 10. To use the package, first install Julia by going to the website: https://julialang.org. Once you have the Julia console open, make sure you have “StatsBase” installed by typing “using Pkg” followed by “Pkg.add("StatsBase”)”. Note there shouldn't be a space between the “add” and the first parenthesis. This is a package required for random number generation. 
Next, download all the files from Github into one folder. At the Julia console type: “include("pathto\\Ground_truth_generation.jl”)” to run the main scrip. Note that “pathto” is the path where you saved the network generation scripts to, and the double backslash is to avoid the Julia misinterpreting the backslash as a control character. A random signaling network will be generated in SBML format called “sampleNetwork.xml”. There are some configuration settings in the Julia script file “Ground_truth_generation.jl”. These include:

1) The number of species “nSpeces” and maximum number of reactions “nRxns_limitation” are modifiable.
2) Randomly assigned ranges for species concentrations and rate constants can be modified via “rnd_species” and “rnd_parameter”.
3) “concentration_perturb” can be used to set the factor that perturbs the concentration at the input species. Default is set to “2”.
4) The number of random networks to generate can be set by changing the variable “sampleSize”.

The file “roadrunner_c_api.dll” is the dynamic library for libRoadRunner that is used to provide SBML simulation support.   Ground_truth_generation.jl refers to https://github.com/Lukez-pi/NetworkGenerator.jl, and rr_funcs-Jin.jl refers to https://github.com/SunnyXu/RoadRunner.jl/tree/master/src. 

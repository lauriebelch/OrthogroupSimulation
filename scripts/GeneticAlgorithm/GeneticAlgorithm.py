# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 11:19:47 2026

@author: JH
"""
import sys
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
import numpy
import pygad
import random
import shutil
import subprocess
import time

import scripts.GeneticAlgorithm.ScoreSomeParameters as ScoreSomeParameters




def GA_workflow(complete_parameters, output_abolsute_path, threads, current_file_path,default_species_pick):
    ###### import scripts here...

    os.mkdir(os.path.join(complete_parameters["output"],"Simulation_Temp_Files"))
    def fitness_func(ga_instance, solution, solution_idx):

        new_solution = []
        for i in solution:
            if i < 0:
                new_solution.append(0.01)
            else:
                new_solution.append(i)
                
                
                
                
                
        print("input = " + str([new_solution[0],new_solution[1],new_solution[2]]))
        try:
            output, outputdir = ScoreSomeParameters.main(new_solution[0],new_solution[1],new_solution[2],complete_parameters,current_file_path,default_species_pick)
        except subprocess.CalledProcessError as e:
                outputdir = os.path.join("/local/home/biol0216/Simulation/Tuning/GA/GA/sim_training_v1/simulation_runs/", [i for i in (str(e).split("/")) if len(i) == 32][0])
                output = "1000,1000,1000"
       	shutil.rmtree(outputdir) 
       	print([new_solution[0],new_solution[1],new_solution[2]])
       	fitness_results = output.split(",")
       	fitness_values = []
       	fitness_results = output.split(",")
       	for i in fitness_results:
            I = 1 + float(i)
            fitness_values.append(1/float(I))
       	print("Results = " + str(fitness_values))
       	return fitness_values
	

    #### read function...

    """
    {'output': 'GA_TEST_OUTPUT', 
     'Orthogroups': '50', 
     'generations': '5', 
     'sample_size': '100', 
     'num_parents_mating': '4', 
     'sol_per_pop': '16', 
     'num_genes': '4', 
     'parent_selection_type': 'nsga2', 
     'keep_parents': '1', 
     'crossover_type': 'uniform', 
     'crossover_probability': '0.9', 
     'mutation_type': 'random', 
     'mutation_percent_genes': '10', 
     'save_solutions': '1', 
     'parallel_processing': '16', 
     'Orthogroup_train_results': 'PARAMETERS_FOLDER'}
    
    """
    
    ## deafult..
    function_inputs = [0,0,0] 

    fitness_function = fitness_func
    num_generations = int(complete_parameters["generations"])
    num_parents_mating = int(complete_parameters["num_parents_mating"])
    sol_per_pop = int(complete_parameters["sol_per_pop"])
    num_genes = len(function_inputs)
    #parent_selection_type = "sss"
    parent_selection_type = complete_parameters["parent_selection_type"]
    keep_parents = int(complete_parameters["keep_parents"])
    crossover_type = complete_parameters["crossover_type"]
    crossover_probability = float(complete_parameters["crossover_probability"])
    mutation_type = complete_parameters["mutation_type"]
    mutation_percent_genes = int(complete_parameters["mutation_percent_genes"])
    sample_size = int(complete_parameters["sample_size"])
    start = time.time()
    
    
    
    ga_instance = pygad.GA(num_generations=num_generations,
                       sample_size = sample_size,
                       num_parents_mating=num_parents_mating,
                       fitness_func=fitness_function,
                       sol_per_pop=sol_per_pop,
                       num_genes=num_genes,
                       parent_selection_type=parent_selection_type,
                       keep_parents=keep_parents,
                       crossover_type=crossover_type,
                       crossover_probability = crossover_probability,
                       mutation_type=mutation_type,
                       mutation_percent_genes=mutation_percent_genes,
                       save_solutions=True,
		       parallel_processing=int(complete_parameters['parallel_processing'])
			)


    end = time.time()
    ga_instance.run()

    solution, solution_fitness, solution_idx = ga_instance.best_solution()
    print("Parameters of the best solution : {solution}".format(solution=solution))
    print("Fitness value of the best solution = {solution_fitness}".format(solution_fitness=solution_fitness))
    print("took : %s" % str(end- start))
    quit()
	#ga_instance.plot_new_solution_rate(save_dir="new_solution_rate")
	#ga_instance.plot_fitness(save_dir="plot_fitness")



















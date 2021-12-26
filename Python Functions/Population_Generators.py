import os
import subprocess
import shutil
import numpy as np
from pathlib import Path
import datetime
import pandas as pd
import csv
import copy
import sys


def generate_population(pop_sze, bounds, constraint):
    population = np.zeros((pop_sze, len(bounds)))
    # generating the initial population
    ii = 0
    while ii < pop_sze:
        for jj in range(len(bounds[:, 0])):
            population[ii, jj] = np.random.uniform(low=bounds[jj, 0], high=bounds[jj, 1], size=1)
        if (constraint and constraint_function(population[ii, 0], population[ii, 1], population[ii, 2]))\
                or not constraint:
            ii += 1

    return population


def import_init_population(init_pop_loc):
    # the imported population must be built a following:
    # first row: [empty]
    # first row: [generation number, convergence]
    # second row: [current winner data]
    # third to last rows: [data of rest of the population]
    init_population = pd.read_excel(init_pop_loc)
    init_population = init_population.to_numpy()
    init_gen = int(init_population[0, 0])
    winner = init_population[1, :]
    init_population = init_population[2:, :]
    return init_population, init_gen, winner


def constraint_function(FW, alpha, beta):
    # alpha - the FD/FW ratio
    # beta the ST/FW ratio
    if FW*((beta ** 2 + alpha ** 2) ** (1/2) - 1) > float(100):
        return True
    else:
        return False


def generate_genetic_population(breeding_pool, mutation_rate, mutation_range, elitism, bounds, constraint):
    elitism = int(elitism)
    par_num = len(breeding_pool[0]) - 1  # the number of optimization parameters
    new_generation = np.zeros((len(breeding_pool), par_num + 1))
    comb_matrix = get_binary_combinations(par_num)
    parents = np.zeros((2, par_num))
    breeding_pool = breeding_pool[breeding_pool[:, par_num].argsort()]
    prob_vec = (1 / breeding_pool[:, par_num]) / sum(1 / breeding_pool[:, par_num])

    ii = 0
    while ii < len(breeding_pool):
        if ii < elitism:  # passing on the fittest among the initial population defined by elitism
            new_generation[ii, :] = breeding_pool[ii, :]
        else:  # the rest of the population is chosen according to a roulette method
            parents[0, :] = breeding_pool[np.argmax(np.random.multinomial(1, prob_vec)), :par_num]
            parents[1, :] = breeding_pool[np.argmax(np.random.multinomial(1, prob_vec)), :par_num]
            jj = np.random.randint(0, len(comb_matrix[:, 0]))
            for kk in range(0, par_num):
                new_generation[ii, kk] = parents[comb_matrix[jj, kk], kk]

            #  inserting mutation according to the mutation rate and range
            for jj in range(0, par_num):
                if np.random.rand() < mutation_rate:
                    bot_bound = max((1 - mutation_range) * new_generation[ii, jj], bounds[jj, 0])
                    top_bound = min((1 + mutation_range) * new_generation[ii, jj], bounds[jj, 1])
                    new_generation[ii, jj] = np.random.uniform(low=bot_bound, high=top_bound)

        if (constraint and constraint_function(new_generation[ii, 0], new_generation[ii, 1], new_generation[ii, 2]))\
                or ~constraint:
            ii += 1

    return new_generation


def get_binary_combinations(param_num):
    comb_mat = np.zeros((2**param_num, param_num))
    for ii in range(2**param_num):
        num = ii
        jj = 0
        while num > 0:
            comb_mat[ii, jj] = num % 2
            num = int(num/2)
            jj += 1

    return comb_mat[1: len(comb_mat[:, 0])-1, :].astype(int)

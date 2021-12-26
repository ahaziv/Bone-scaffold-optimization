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
    # population[0, :] = [0.563, 2, 0.8]        # in case you want to evaluate a certain geometry

    population = np.zeros((pop_sze, len(bounds)))
    # generating the initial population
    ii = 0
    while ii < pop_sze:
        for jj in range(len(bounds[:, 0])):
            population[ii, jj] = np.random.uniform(low=bounds[jj, 0], high=bounds[jj, 1], size=1)

        if (constraint and constraint_function(population[ii, 0], population[ii, 1], population[ii, 2]))\
                or ~constraint:
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


def create_mesh(FW, FD, ST, ES):
    # writing a command script for ANSYS creating the mesh
    create_mesh_path = os.path.join(os.getcwd(), 'Create_Mesh.wbjn')
    # cleaning the mesh Mesh_files directory
    clean_file = os.path.join(os.getcwd(), 'Mesh_files\dp0\SYS\MECH\SYS.dat')
    if os.path.isfile(clean_file):
        os.remove(clean_file)

    with open(create_mesh_path, 'w') as f:
        f.write('SetScriptVersion(Version="19.0.136")\n'
                'Open(FilePath=GetAbsoluteUserPathName("' + os.getcwd().replace("\\", "/") + '\Mesh.wbpj"))\n'
                'designPoint1 = Parameters.GetDesignPoint(Name="0")\n'
                'parameter1 = Parameters.GetParameter(Name="P1")   # Fillament Width\n'
                'designPoint1.SetParameterExpression(\n'
                '    Parameter=parameter1,\n'
                '    Expression="' + str(FW) + ' [mm]")\n'
                'parameter2 = Parameters.GetParameter(Name="P2")   # Fillament Distance\n'
                'designPoint1.SetParameterExpression(\n'
                '    Parameter=parameter2,\n'
                '    Expression="' + str(FD) + ' [mm]")\n'
                'parameter3 = Parameters.GetParameter(Name="P3")   # Slice Thickness\n'
                'designPoint1.SetParameterExpression(\n'
                '    Parameter=parameter3,\n'
                '    Expression="' + str(ST) + ' [mm]")\n'
                'parameter4 = Parameters.GetParameter(Name="P4")   # Element Size\n'
                'designPoint1.SetParameterExpression(\n'
                '    Parameter=parameter4,\n'
                '    Expression="' + str(ES) + ' [m]")\n'
                # 'Save(Overwrite=True)\n'
                'Update()\n')
        f.close()

    # running ANSYS in batch mode while executing the CreateMeshDir script
    run_ANS_command = '"D:\programs\ANSYS Inc\\v190\Framework\\bin\Win64\\runwb2"' + ' -R ' + '"' + create_mesh_path + '"'
    subprocess.run(run_ANS_command)


def reset_model():
    os.remove(os.path.join(os.getcwd(), 'Mesh.wbpj'))
    shutil.rmtree(os.path.join(os.getcwd(), 'Mesh_files'), ignore_errors=True)
    run_ANS_command = '"D:\programs\ANSYS Inc\\v190\Framework\\bin\Win64\\runwb2"' + ' -R ' + '"' +\
                      os.getcwd() + '\Backup_Model\Reset_Model.wbjn"'

    subprocess.run(run_ANS_command)


def calc_loss(FW, FD, ST):
    # this function imports/creates two input files for a MATLAB code and then runs it.
    # One input file includes FW, FD and ST. The other file is a mesh data file,
    # if the second file is nonexistant or corrupt the function will return the
    # value [0,0,10^7] as an output. otherwise the three loss values will be returned in an array.

    # moving the mesh file to the main directory and deleting it form the former
    source = os.path.join(os.getcwd(), 'Mesh_files\dp0\SYS\MECH\SYS.dat')
    destination = os.path.join(os.getcwd(), 'SYS.dat')

    if os.path.isfile(source):
        if os.path.isfile(destination):
            os.remove(destination)
        shutil.move(source, destination)
    else:
        return [0, 0, 10**6]

    # creating a txt file containing the parameters
    data_file_path = os.path.join(os.getcwd(), 'Parameters.txt')
    data_file = open(data_file_path, "w")
    data_file.write(str(FW) + '\n' + str(FD) + '\n' + str(ST))
    data_file.close()
    # writing the loss data file
    loss_file_path = os.path.join(os.getcwd(), 'Loss.txt')
    # running the matlab script
    run_matlab_command = '"D:\programs\\bin\matlab.exe"' + ' -nosplash -nodesktop -wait -r ' + \
                         '"run(\'' + os.getcwd() + '\Matlab_Script\CCTG_OPT_main_function.m\');exit;"'

    p = subprocess.Popen(run_matlab_command)
    stdout, stderr = p.communicate()
    loss_file = open(loss_file_path, "r")
    curve_loss = loss_file.readline()
    curve_loss.replace('\n', '')
    pore_loss = loss_file.readline()
    pore_loss.replace('\n', '')
    loss = loss_file.readline()
    loss.replace('\n', '')
    loss_arr = [float(curve_loss), float(pore_loss), float(loss)]
    return loss_arr


def write_data(file_location, file_name, individual, is_winner):
    if individual[5] != 10 ** 6 and individual[5] != 10 ** 5:
        if is_winner:
            with open(os.path.join(file_location, file_name + '.csv'), 'a', newline='') as csvfile:
                filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                filewriter.writerow(["{:.{}f}".format(individual[0], precision),
                                     "{:.{}f}".format(individual[1], precision),
                                     "{:.{}f}".format(individual[2], precision),
                                     "{:.{}f}".format(individual[3], precision + 2),
                                     "{:.{}f}".format(individual[4], precision + 2),
                                     "{:.{}f}".format(individual[5], precision + 2), '1'])

        else:
            with open(os.path.join(file_location, file_name + '.csv'), 'a', newline='') as csvfile:
                filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
                filewriter.writerow(["{:.{}f}".format(individual[0], precision),
                                     "{:.{}f}".format(individual[1], precision),
                                     "{:.{}f}".format(individual[2], precision),
                                     "{:.{}f}".format(individual[3], precision + 2),
                                     "{:.{}f}".format(individual[4], precision + 2),
                                     "{:.{}f}".format(individual[5], precision + 2), '0'])

    else:
        if individual[5] == 10 ** 6:
            failure_type = 'Geometrical Failure'
        else:
            failure_type = 'Algorithm Failure'
        with open(os.path.join(file_location, file_name + '.csv'), 'a', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(["{:.{}f}".format(individual[0], precision),
                                 "{:.{}f}".format(individual[1], precision),
                                 "{:.{}f}".format(individual[2], precision), failure_type])


def write_gen_data(file_location, file_name, gen_num, algo_flag, convergence, bounds):
    if algo_flag == 'greedy':
        with open(os.path.join(file_location, file_name + '.csv'), 'a', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(['Generation:', str(gen_num), 'Convergance rate:', "{:.{}f}".format(convergence, precision)])
            filewriter.writerow(["{:.{}f}".format(bounds[0, 0], precision), "{:.{}f}".format(bounds[0, 1], precision),
                                 "{:.{}f}".format(bounds[1, 0], precision), "{:.{}f}".format(bounds[1, 1], precision),
                                 "{:.{}f}".format(bounds[2, 0], precision), "{:.{}f}".format(bounds[2, 1], precision), '2'])

    if algo_flag == 'genetic':
        with open(os.path.join(file_location, file_name + '.csv'), 'a', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(['Generation:', str(gen_num)])


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



##################################### main block #####################################################
flag_algorithm = 'genetic'  # can be either 'greedy' or 'genetic'
flag_import_pop = 1

ES_factor = 5 * 10 ** -5    # a factor to determine the average element size
precision = 8               # the number of digits printed out
par_num = 3                 # number of optimization parameters
pop_sze = 5
gen_num = 196                # the number of generations to run

bounds = np.array([[55, 750], [1.1, 3], [0.6, 0.9]])
search_bnds = copy.deepcopy(bounds)

flag_test = False
constraint = True
# parameter boundaries
# The first row of bounds represents the FW boundaries in [mm]
# The second row of bounds represents FD/FW boundaries [unitless]
# The Third row of bounds represents ST/FW boundaries [unitless]

if flag_algorithm == 'greedy':
    init_pop_sze = 8
    convergence = 1
    conv_rate = 0.978                # the convergence rate at each generation (greedy algorithm)
    winner = np.mean(bounds, axis=1)
    winner = np.append(winner, [0, 0, 100])

elif flag_algorithm == 'genetic':
    mutation_rate = 0.546        # the mutation rate of an individual at each generation (genetic algorithm)
    mutation_range = 0.1
    elitism = 4
    init_pop_sze = pop_sze
    conv_rate = 1

else:
    print('Invalid algorithm type.')
    sys.exit()

population = np.zeros((gen_num, init_pop_sze, par_num + 3))

# generating/importing the initial population
if flag_import_pop:
    imp_pop_loc = os.path.join(os.getcwd(), 'Data Files\initial_population.xlsx')
    imported_pop, init_gen_num, winner = import_init_population(imp_pop_loc)
    if init_gen_num != 0:
        init_pop_sze = copy.deepcopy(pop_sze)

    imp_pop_sze = len(imported_pop[:, 0])
    if (flag_algorithm == 'genetic') and imp_pop_sze != init_pop_sze:
        raise ValueError('Can\'t import a partial generation into a genetic algorithm')

    population[init_gen_num, :imp_pop_sze, :] = imported_pop
    if flag_algorithm == 'greedy':
        convergence = convergence * conv_rate ** init_gen_num
        search_bnds_int = (bounds[:, 1] - bounds[:, 0]) * convergence
        for kk in range(3):
            search_bnds[kk, 0] = max(bounds[kk, 0], winner[kk] - convergence * search_bnds_int[kk] / 2)
            search_bnds[kk, 1] = min(bounds[kk, 1], winner[kk] + convergence * search_bnds_int[kk] / 2)
        if imp_pop_sze < init_pop_sze:
            population[init_gen_num, imp_pop_sze:init_pop_sze, :3] = generate_population(init_pop_sze - imp_pop_sze, search_bnds)

    # if flag_algorithm == 'genetic':
    #     population[init_gen_num + 1, :, :] =\
    #         generate_genetic_population(imported_pop, mutation_rate, mutation_range, elitism, bounds)
else:
    imp_pop_sze = 0
    init_gen_num = 0
    population[0, :, 0:3] = generate_population(init_pop_sze, bounds, constraint)

# writing 2 files, one to store data and the other to record failed runs,
# file names are generated from the date and hour/minute of the run
file_location = os.path.join(os.getcwd(), 'Data Files')
data_file_name = 'Data' + datetime.datetime.now().strftime("%d-%m-%Y %H-%M")
with open(os.path.join(file_location, data_file_name + '.csv'), 'a', newline='') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(['FW', 'FD/FW', 'ST/FW', 'curve_loss', 'pore_loss', 'tot_loss', 'conv_rate',
                         str(conv_rate)])
if flag_algorithm == 'greedy':
    write_gen_data(file_location, data_file_name, init_gen_num, flag_algorithm, convergence, search_bnds)
    write_data(file_location, data_file_name, winner, 1)
else:
    write_gen_data(file_location, data_file_name, init_gen_num, flag_algorithm, 0, bounds)

failure_file_name = 'Failure Data' + datetime.datetime.now().strftime("%d-%m-%Y %H-%M")
with open(os.path.join(file_location, failure_file_name + '.csv'), 'a', newline='') as csvfile:
    filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
    filewriter.writerow(['FW', 'FD/FW', 'ST/FW', 'Failure_Type'])

# running the optimization process
for ii in range(init_gen_num, gen_num):
    if flag_algorithm == 'greedy' and ii == 0:
        jj_max = init_pop_sze
    else:
        jj_max = pop_sze
    jj = 0
    if imp_pop_sze:
        jj = imp_pop_sze
        imp_pop_sze = 0
    while jj < jj_max:
        if ii != 0 and flag_algorithm == 'genetic' and jj < elitism:
            population[ii, jj, 5] = population[ii, jj, 3]
            population[ii, jj, 3] = 0
        else:
            # the individual properties
            FW = population[ii, jj, 0]  # Filament Width
            FD = population[ii, jj, 1]*population[ii, jj, 0]  # Filament Distance
            ST = population[ii, jj, 2]*population[ii, jj, 0]  # Slice Thickness
            ES = (FW + FD + FW) * ES_factor/3  # Element Size
            # a test to see if the algorithm converges to: (FW = 45, FD/FW = 2, ST/FW = 0.8)
            if flag_test:
                population[ii, jj, 5] = ((45 - FW) ** 2 + (45 / 2 * (2 - FD/FW)) ** 2 + (
                        45 / 0.8 * (0.8 - ST/FW)) ** 2) ** 0.5
            else:
                # creating the SYS.dat mesh file for the CCTG algorithm
                create_mesh(FW, FD, ST, ES)
                # running the CCTG algorithm to calculate the appropriate loss
                population[ii, jj, 3:] = calc_loss(FW, FD, ST)
        if population[ii, jj, 5] != 10 ** 6 and population[ii, jj, 5] != 10 ** 5:
            write_data(file_location, data_file_name, population[ii, jj, :], 0)
            jj += 1
        else:
            write_data(file_location, failure_file_name, population[ii, jj, :], 0)
            if population[ii, jj, 5] == 10**6:
                reset_model()
            if flag_algorithm == 'greedy' or ii == 0:
                population[ii, jj, 0:3] = generate_population(1, search_bnds, constraint)
            else:
                population[ii, jj, 0:5] = population[ii-1, jj, 0:5]
                continue

    if ii < gen_num-1:
        if flag_algorithm == 'greedy':
            convergence = convergence*conv_rate
            min_index = np.argmin(population[ii, :jj_max, 5])
            if population[ii, min_index, 5] < winner[5]:
                winner = population[ii, min_index, :]

            search_bnds_int = (bounds[:, 1] - bounds[:, 0])*convergence
            for kk in range(3):
                search_bnds[kk, 0] = max(bounds[kk, 0], winner[kk] - search_bnds_int[kk]/2)
                search_bnds[kk, 1] = min(bounds[kk, 1], winner[kk] + search_bnds_int[kk]/2)

            population[ii + 1, :pop_sze, :3] = generate_population(pop_sze, search_bnds, constraint)
            write_gen_data(file_location, data_file_name, ii + 1, flag_algorithm, conv_rate ** (ii + 1), search_bnds)
            write_data(file_location, data_file_name, winner, 1)

        elif flag_algorithm == 'genetic':
            send_dat = np.delete(population[ii, :, :], [3, 4], axis=1)
            population[ii + 1, :, :4] =\
                generate_genetic_population(send_dat, mutation_rate, mutation_range, elitism, bounds, constraint)
            write_gen_data(file_location, data_file_name, ii + 1, flag_algorithm, 0, bounds)

    print('Generation Number: ' + str(ii) + '\nAverage Loss: ' + str(np.mean(population[ii, :jj_max, 5])))

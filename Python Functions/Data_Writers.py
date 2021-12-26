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


def write_data(file_location, file_name, individual, is_winner, precision):
    if individual[len(individual)-1] != 10 ** 6 and individual[len(individual)-1] != 10 ** 5:
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
            failure_type = 'CCTG Algorithm Failure'
        with open(os.path.join(file_location, file_name + '.csv'), 'a', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(["{:.{}f}".format(individual[0], precision),
                                 "{:.{}f}".format(individual[1], precision),
                                 "{:.{}f}".format(individual[2], precision), failure_type])


def write_gen_data(file_location, file_name, gen_num, algo_flag, convergence, bounds, precision):
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

    if algo_flag == 'pattern':
        with open(os.path.join(file_location, file_name + '.csv'), 'a', newline='') as csvfile:
            filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
            filewriter.writerow(['Generation:', str(gen_num)])

def write_bounds(file_location, file_name, bounds, precision):
    with open(os.path.join(file_location, file_name + '.csv'), 'a', newline='') as csvfile:
        filewriter = csv.writer(csvfile, delimiter=',', quotechar='|', quoting=csv.QUOTE_MINIMAL)
        filewriter.writerow(["{:.{}f}".format(bounds[0, 0], precision), "{:.{}f}".format(bounds[0, 1], precision),
                             "{:.{}f}".format(bounds[1, 0], precision), "{:.{}f}".format(bounds[1, 1], precision),
                             "{:.{}f}".format(bounds[2, 0], precision), "{:.{}f}".format(bounds[2, 1], precision), '2'])


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
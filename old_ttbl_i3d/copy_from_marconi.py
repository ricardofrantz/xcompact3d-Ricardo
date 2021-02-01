from os import getcwd, system
from os.path import normpath, basename
import time

workdir = basename(normpath(getcwd()))
path = "lmortari@login.marconi.cineca.it:/marconi_work/IscrC_SUB-SAB_0/" + workdir + "/"
print(path)

file = []
file.append(['', '*.f90', '*.prm', '*.xdmf', 'Makefile','logfile'])
file.append(['/out/', '*'])
file.append(['/stats/', '*'])
file.append(['/data/', 'phi1*'])
file.append(['/data/', '*'])
file.append(['/*/', '*'])

for i in file:
    folder = i.pop(0)
    print(folder)
    for j in i:
        print(j)
        system('rsync -auv '+path+folder+j+' .'+folder)

while(True):

    file = []
    #file.append(['', '*.out'])
    #file.append(['', '*.f90', '*.prm', '*.xdmf', '*.out', 'Makefile'])
    #file.append(['/data/', '*m????', '*m?????', '*????'])
    #file.append(['/out/', '*'])
    #file.append(['/stats/', '*'])

    for i in file:
        folder = i.pop(0)
        print(folder)
        for j in i:
            print(j)
            system('rsync -auv '+path+folder+j+' .'+folder)
    
    time.sleep(450)

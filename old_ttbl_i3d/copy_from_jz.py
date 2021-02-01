from os import getcwd, system
from os.path import normpath, basename
import time

workdir = basename(normpath(getcwd()))
path = "rvpo014@jean-zay.idris.fr:/gpfswork/rech/vpo/rvpo014/side/ttbl/ttbl400"
print(path)

file = []
file.append(['/data/','utmap*'])
file.append(['/out/', '*'])
file.append(['/', 'Makefile','*.prm','*.f90', '*.png', '*.py', 'logfile','s_*','theta','yp.dat','*.xdmf'])

for i in file:
    folder = i.pop(0)
    print(folder)
    for j in i:
        print(j)
        system('rsync -auv '+path+folder+j+' .'+folder)

import numpy as np
import subprocess
import matplotlib.pyplot as plt
import pdb
import os

# Parameters
# TODO adapt to what you need (folder path executable input filename)
executable = 'Exercice1_student'  # Name of the executable (NB: .exe extension is required on Windows)
repertoire = r"/Users/Sayu/Desktop/3-body_problem_physNum/"
os.chdir(repertoire)

input_filename = 'configuration.in.example'  # Name of the input file

alpha = np.array([0,0.5,1])
nsteps = np.array([4000, 6000, 10000, 14e3, 20e3])
nsimul = len(alpha)  # Number of simulations to perform

tfin = 259200  # TODO: Verify that the value of tfin is EXACTLY the same as in the input file

dt = tfin / nsteps

paramstr = 'alpha'  # Parameter name to scan
param = alpha  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
convergence_list = []

for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    cmd = f"./{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    cmd = f"./{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

error = np.zeros(nsimul)

fig, ax = plt.subplots(constrained_layout=True)
plt.grid(True, linestyle="--", alpha=0.3)
fig2, ax2 = plt.subplots(constrained_layout=True)
plt.grid(True, linestyle="--", alpha=0.3)
fig3, ax3 = plt.subplots(constrained_layout=True)
plt.grid(True, linestyle="--", alpha=0.3)
for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    vx = data[-1, 1]  # final position, velocity, energy
    vy = data[-1, 2]
    xx = data[-1, 3]
    yy = data[-1, 4]
    En = data[-1, 5]
    convergence_list.append(xx)
    # TODO compute the error for each simulation

    x0=data[0,3]
    y0=data[0,4]

    error[i] =  np.abs(yy)
    lw = 1.5
    fs = 16

    
    ax.plot(data[:, 3], data[:, 4])
    ax2.plot(data[:, 3], data[:, 4])
    #fig3, ax3 = plt.subplots(constrained_layout=True)
    plt.grid(True, linestyle="--", alpha=0.3)
    ax3.plot(t,data[:,5])
    ax3.set_xlabel('$t$ [s]', fontsize=fs)
    ax3.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax3.xaxis.set_major_formatter(plt.ScalarFormatter(useMathText=True))
    # Appliquer la notation scientifique à l'axe X
    ax3.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
    ax3.set_ylabel('$E_{mec}$ [J]', fontsize=fs)

# Ajouter un disque
disque = plt.Circle((3.80321e+08, 0), 1737100, color='black')
disque2 = plt.Circle((3.80321e+08, 0), 1737100, color='black')
disque3 = plt.Circle((341878931, 0), 500000, color='red')
# Ajouter le disque au graphe
ax.add_patch(disque)
ax.add_patch(disque3)
ax2.add_patch(disque2)
ax.set_xlabel('$x\'$ [m]', fontsize=fs)
ax.set_ylabel('$y\'$ [m]', fontsize=fs)
ax2.set_xlabel('$x\'$ [m]', fontsize=fs)
ax2.set_ylabel('$y\'$ [m]', fontsize=fs)

# Ajuster les limites du graphique
ax.set_xlim(3.4e8,3.9e8)
ax.set_ylim(-1e7,0.31e8)
ax2.set_xlim(3.7e8,3.85e8)
ax2.set_ylim(-6e6,6e6)
ax.set_aspect('equal')  # Assurer un aspect circulaire
ax2.set_aspect('equal')  



#    error[i] =  np.sqrt((yy-y0)**2+(xx-x0)**2) # quoi prendre ?? 

# uncomment the following if you want debug
#import pdb
#pbd.set_trace()
'''
plt.figure()
plt.loglog(dt, error, 'r+-', linewidth=lw)
plt.xlabel('dt [s]', fontsize=fs)
plt.ylabel('final position error [m]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)



a = (error[-1]-error[-2])/(dt[-1]-dt[-2])

b = error[-1] - a*dt[-1]

error -= b

plt.figure()
plt.loglog(dt, -error, 'r+-', linewidth=lw)
plt.xlabel('dt [s]', fontsize=fs)
plt.ylabel('final position error [m]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

"""
Si on n'a pas la solution analytique: on représente la quantite voulue
(ci-dessous v_y, modifier selon vos besoins)
en fonction de (Delta t)^norder, ou norder est un entier.
"""
norder = 1  # Modify if needed

a = (convergence_list[-1]-convergence_list[-2])/((dt[-1]**norder)-(dt[-2]**norder))

b = convergence_list[-1] - a*(dt[-1]**norder)

plt.figure()
plt.plot(dt**norder, convergence_list, 'k+-', linewidth=lw)
plt.xlabel('dt [s]', fontsize=fs)
plt.ylabel('xx [m]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

plt.figure()
plt.loglog(dt**norder, convergence_list-b, 'r+-', linewidth=lw)
plt.xlabel('dt [s]', fontsize=fs)
plt.ylabel('final position error [m]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)

#plt.figure()
#plt.plot(t, data[:,5], 'k+-', linewidth=lw)
#plt.xlabel('t [s]', fontsize=fs)
#plt.ylabel('Énergie mécanique [J]', fontsize=fs)
#plt.xticks(fontsize=fs)
#plt.yticks(fontsize=fs)
#plt.grid(True)


'''
print(nsteps)
plt.show()

#%% Libraries

import numpy as np
import matplotlib.pyplot as plt
import math

#%% Parameters
k = 10 #erg/cm2
x0 = 3  #cm, it is the equilibrium length of the spring 
m1 = 10 #g
m2 = 12 #g
simulation_time = 15

#%% Initial Conditions

x1_0 = 2 #cm
x2_0 = 6 #cm
v1_0 = 0 # cm/s
v2_0 = 0 #cm /s 
F_initial = k *(x2_0 - x1_0 - x0) # g * (cm /s^2)

# %% Velocity Verlet algorithm

def VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial):
    
    time = np.arange(0, simulation_time + delta_t, delta_t)
    
    x1 = [x1_0]
    x2 = [x2_0]
    v1 = [v1_0]
    v2 = [-v2_0] #the negative sign is because the direction is opposite
    
    F = [F_initial]

    U = [(1/2) * k * (x2[0] - x1[0]- x0)**2]
    energy_kinetic_1 =[(1/2) * m1 * (v1_0)**2]
    energy_kinetic_2 = [(1/2) * m2 * (v2_0)**2]

    total_energy = [(U[0] + energy_kinetic_1[0] + energy_kinetic_2[0])]

    i = 0
    for i in range(len(time)-1):
            x1.append( x1[i] + delta_t * v1[i] + delta_t**2 * F[i]/(2*m1) )
            x2.append( x2[i]+ delta_t * v2[i] - delta_t**2 * F[i]/(2*m2) )

            F.append( k * (x2[i+1] - x1[i+1]- x0) )

            v1.append( v1[i] + (delta_t * (F[i] + F[i+1])) / (2*m1) )
            v2.append( v2[i] - (delta_t * (F[i] + F[i+1])) / (2*m2) )

            U.append((1/2) * k * (x2[i+1] - x1[i+1]- x0)**2) 
            energy_kinetic_1.append((1/2) * m1 * (v1[i+1])**2)
            energy_kinetic_2.append((1/2) * m2 * (v2[i+1])**2)
            
            total_energy.append(U[i+1]+ energy_kinetic_1[i+1] + energy_kinetic_2[i+1])

            i += 1

    return x1, x2, v1, v2, F, time, U, energy_kinetic_1, energy_kinetic_2, total_energy
           #0  #1  #2  #3  #4  #5  #6       #7                #8               #9

#%% 3)
delta_t_test = [10**-5, 10**-4, 10**-3, 10**-2, 10**-1, 1]

for i in range(len(delta_t_test)-1):
        total_energy = VelocityVerletAlgorithm(simulation_time, delta_t_test[i], x0, x1_0, x2_0, v1_0, v2_0, F_initial)[9]
        time = VelocityVerletAlgorithm(simulation_time, delta_t_test[i], x0, x1_0, x2_0, v1_0, v2_0, F_initial)[5]
        plt.ylim(4.7, 5.1) #to set-up depend your result
        plt.plot(time, total_energy, '--',label=" Δt = {}".format(delta_t_test[i]))
        plt.legend(fontsize=13)
        plt.xlabel('Time')
        plt.ylabel('Total Energy')
        plt.title("HW1 - Point 3 - Total Energy at different Δt")      
plt.show()

# %% 4) Make plots of the particle coordinates and velocities versus time based on the simulation results.

delta_t = 0.01

x1_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[0]
x2_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[1]
v1_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[2]
v2_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[3]
time = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[5]

plt.plot(time, x1_verlet, 'o', color ='b',label="x1")
plt.plot(time, x2_verlet, 'o', color ='r',label="x2")
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Particle cordinates')
plt.title("HW1 - Point 4 - Velocity Verlet Algorithm")
plt.show()

plt.plot(time, v1_verlet, 'o', color ='b',label="v1")
plt.plot(time, v2_verlet, 'o', color ='r',label="v2")
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Velocity ')
plt.title("HW1 - Point 4 - Velocity Verlet Algorithm")
plt.show()

# %% 5)
total_energy_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[9]
plt.plot(time, total_energy_verlet , '-' ,color ='red', label="Total energy")
plt.ylim(4, 6)
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title("HW1 - Point 4 - Velocity Verlet Algorithm")
plt.show()

U_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[6]
plt.plot(time, U_verlet,'--',label="Potential Energy")
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title("HW1 - Point 4 - Velocity Verlet Algorithm")


energy_kinetic1_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[7]
energy_kinetic2_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[8]
plt.plot(time, energy_kinetic1_verlet,'--', label="Kinetic Energy particule 1")
plt.plot(time, energy_kinetic2_verlet,'--',label="Kinetic Energy particule 2")
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title("HW1 - Point 4 - Velocity Verlet Algorithm")
plt.show()

#%% 6)

#For time steps greater than 1 the vibrational frequency is not uniform and decrease, 
# moreover, the plot shape starts to have distortions and is not uniform.
# For small stiffness of the spring the frequency decrease considerably compared with high values,
# additionally for bigger stiffness of the spring, the amplitude of the wave increase through time. 

#%% 7)
velocity_light  = 29979245800  #cm/s
reduced_mass = (m1 + m2)/2 #g
vibrational_frecuency  = 1/(velocity_light*2*math.pi) * math.sqrt(k/reduced_mass)
print('Vibrational frecuency = ',vibrational_frecuency,  'cm ^-1')

#The analitycal result is simillar with the frecuency of the figures.

#%% 8)

def Euler_Algorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial):
    
    time = np.arange(0, simulation_time + delta_t, delta_t)
    
    x1 = [x1_0]
    x2 = [x2_0]
    v1 = [v1_0]
    v2 = [-v2_0]
    
    F = [F_initial]

    U = [(1/2) * k * (x2[0] - x1[0]- x0)**2]
    energy_kinetic_1 =[(1/2) * m1 * (v1_0)**2]
    energy_kinetic_2 = [(1/2) * m2 * (v2_0)**2]

    total_energy = [(U[0] + energy_kinetic_1[0] + energy_kinetic_2[0])]

    i = 0
    for i in range(len(time)-1):
            x1.append( x1[i] + delta_t * v1[i] + delta_t**2 * F[i]/(2*m1) )
            x2.append( x2[i]+ delta_t * v2[i] - delta_t**2 * F[i]/(2*m2) )

            F.append( k * (x2[i+1] - x1[i+1]- x0) )

            v1.append( v1[i] + (delta_t * (F[i] / (m1) )))
            v2.append( v2[i] - (delta_t * (F[i] / (m2) )))

            U.append((1/2) * k * (x2[i+1] - x1[i+1]- x0)**2) 
            energy_kinetic_1.append((1/2) * m1 * (v1[i+1])**2)
            energy_kinetic_2.append((1/2) * m2 * (v2[i+1])**2)
            
            total_energy.append(U[i+1]+ energy_kinetic_1[i+1] + energy_kinetic_2[i+1])

            i += 1

    return x1, x2, v1, v2, F, time, U, energy_kinetic_1, energy_kinetic_2, total_energy
           #0  #1  #2  #3  #4  #5  #6       #7                #8               #9

x1_euler = Euler_Algorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[0]
x2_euler = Euler_Algorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[1]

plt.plot(time, x1_euler, 'o', color ='b',label="x1")
plt.plot(time, x2_euler, 'o', color ='r',label="x2")
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Particle cordinates')
plt.title("HW1 - Point 8 - Euler Algorithm")
plt.show()

total_energy_euler = Euler_Algorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[9]
plt.ylim(0, 7)
plt.plot(time, total_energy_euler, '-' ,color ='red',label="Total energy")
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title("HW1 - Point 8 - Euler Algorithm")
plt.show()

#%% Comparation 

plt.plot(time, x1_verlet, '--', color ='b',label="x1-Velocity Verlet Algorithm")
plt.plot(time, x2_verlet, '--', color ='b',label="x2-Velocity Verlet Algorithm")
plt.plot(time, x1_euler, '--', color ='r',label="x1-Euler Algorithm")
plt.plot(time, x2_euler, '--', color ='r',label="x2-Euler Algorithm")
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Particle cordinates')
plt.title("HW1 - Point 8 - Comparison Trayectories")
plt.show()

total_energy_verlet = VelocityVerletAlgorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[9]
total_energy_euler = Euler_Algorithm(simulation_time, delta_t, x0, x1_0, x2_0, v1_0, v2_0, F_initial)[9]
plt.plot(time, total_energy_verlet, '-' ,color ='b', label="Total energy-Velocity Verlet Algorithm")
plt.plot(time, total_energy_euler, '-' ,color ='r',label="Total energy-Euler Algorithm")
plt.ylim(4, 7)
plt.legend(fontsize=9)
plt.xlabel('Time')
plt.ylabel('Energy')
plt.title("HW1 - Point 8 - Comparison Energy Conservation")
plt.show()

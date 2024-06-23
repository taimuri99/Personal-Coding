from Ising_obj import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# for subsequent plotting of the spin glass
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(20, 12))

kBT = 1
side = 30
mc_step = 300

lat = SpinGlass(kBT, side)
print(lat.M)
print(lat.E)
fig.clf()
ax = fig.add_subplot(1, 1, 1)
im_content = ax.imshow(lat.config, cmap=cm.afmhot, interpolation='nearest')
fig_name = 'initial'  + '.png'
plt.savefig(fig_name, bbox_inches='tight')

# Thermalization
for i in range(mc_step):
    # metropolis sweep
    lat.metropolis_flipper()
    print("Magnetization = ", lat.M)
    print("Energy = ", lat.E)
    # add the writing onto a file

    # plotting
    fig.clf()
    ax = fig.add_subplot(1, 1, 1)
    im_content = ax.imshow(lat.config, cmap=cm.afmhot, interpolation='nearest')
    fig_name = str(i//1) + '.png'
    if i % 1 == 0:
        plt.savefig(fig_name, bbox_inches='tight')

# Averages at equilibrium
# Add here the calculation of the specific heat
# calculate the average of energy and energy-squared

import matplotlib.pyplot as plt
import numpy as np

x = np.arange(-5, 5, 0.01)
y = np.arange(-5, 5, 0.01)

x, y= np.meshgrid(x, y)

z = x * y - x + y

fig = plt.figure()

ax = plt.axes(projection='3d')

ax.plot_surface(x,y,z,
                cmap='viridis')
plt.show()
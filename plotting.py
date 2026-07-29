import numpy as np
import matplotlib.pyplot as plt

def P_plot(geo, Fields):

    import numpy as np
    import matplotlib.pyplot as plt

    x = np.linspace(0, geo.lx, geo.Nx)
    y = np.linspace(0, geo.ly, geo.Ny)
    X, Y = np.meshgrid(x, y, indexing='ij')

    plt.figure(figsize=(6,5))

    plt.contourf(X, Y, Fields.P,
                levels=30,
                cmap='coolwarm')

    plt.colorbar(label='Pressure (Pa)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Pressure Contours')
    plt.axis('equal')
    plt.show()
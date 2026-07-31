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

    u_center = 0.5*(Fields.u[1:,:] + Fields.u[:-1,:])
    v_center = 0.5*(Fields.v[:,1:] + Fields.v[:,:-1])

    xv = np.linspace(0, geo.lx, geo.Nx)
    yv = np.linspace(0, geo.ly, geo.Ny)

    Xv, Yv = np.meshgrid(xv, yv, indexing='ij')

    step = 3

    plt.quiver(Xv[::step, ::step],
                Yv[::step, ::step],
                u_center[::step, ::step],
                v_center[::step, ::step],
                scale=0.1,
                color='black')

    plt.grid 
    plt.colorbar(label='Pressure (Pa)')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.title('Pressure Contours')
    plt.axis('equal')
    plt.show()
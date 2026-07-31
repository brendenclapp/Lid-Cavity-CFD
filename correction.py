import numpy as np


def correction(geo, Faces, Fields):

 
    #Pressure
    for i in range (geo.Nx):
        for j in range (geo.Ny):

            Fields.P[i,j] = Fields.P[i,j] + 0.2*Fields.Pp[i,j] 

  
    #X-momentum
    for i in range (1, geo.Nx):
        for j in range (geo.Ny): 

            Fields.u[i,j] = Fields.u[i,j] + 0.8*(Fields.Pp[i-1,j]-Fields.Pp[i,j])*(geo.dy/Faces.a_P_u[i,j]) 

    Fields.u[-1, :] = Fields.u[-2, :]

    
    #Y-momentum
    for i in range (geo.Nx):
        for j in range(1, geo.Ny):

            Fields.v[i,j] = Fields.v[i,j] + 0.8*(Fields.Pp[i,j-1]-Fields.Pp[i,j]) * (geo.dx/Faces.a_P_v[i,j])



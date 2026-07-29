import numpy as np


def correction(geo, var, Faces, Fields, Coupler, Tri):

 
    #Pressure
    for i in range (geo.Nx):
        for j in range (geo.Ny):

            Fields.P[i,j] = Fields.P[i,j] + 0.2*Fields.Pp[i,j] 

    i = 0
    j = 0

    for i in range (1, geo.Nx):
        for j in range (geo.Ny): 

            Fields.u[i,j] = Fields.u[i,j] + 0.8*(Fields.Pp[i-1,j]-Fields.Pp[i,j])*(geo.dy/Faces.a_P_u[i,j]) 

    Fields.u[-1, :] = Fields.u[-2, :]

    i = 0
    j = 0

    for i in range (geo.Nx):
        for j in range(1, geo.Ny):

            Fields.v[i,j] = Fields.v[i,j] + 0.8*(Fields.Pp[i,j-1]-Fields.Pp[i,j]) * (geo.dx/Faces.a_P_v[i,j])



    fmt = {'float_kind': lambda x: f"{x:12.6f}"}

    #print('v final results')
    #print(np.array2string(np.flipud(Fields.v.T),formatter=fmt,max_line_width=100000000000))

   # print('u final results')
    #print(np.array2string(np.flipud(Fields.u.T),formatter=fmt,max_line_width=100000000000))

   # print('P final results')
   # print(np.array2string( np.flipud(Fields.P.T), formatter=fmt, max_line_width=100000000000 ))
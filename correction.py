import numpy as np


def correction(geo, var, Faces, Fields, Coupler, Tri):

    i = 0
    j = 0
    print('P_prime results')
    print(np.array2string(np.flipud(Fields.P_prime.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))

    for i in range (geo.Nx):
        for j in range (geo.Ny):

            Fields.P[i,j] = Fields.P[i,j] + 0.3*Fields.P_prime[i,j] 

    i = 0
    j = 0

    for i in range (1, geo.Nx):
        for j in range (geo.Ny):

            Fields.u_psu[i,j] = 0.7*Fields.u_psu[i,j] + (1 - 0.7)*Fields.u_old[i,j]   

            Fields.u[i,j] = Fields.u_psu[i,j] + (Coupler.d_we[i,j]*(Fields.P_prime[i-1,j]-Fields.P_prime[i,j])) 

    Fields.u[-1, :] = Fields.u[-2, :]

    i = 0
    j = 0

    for i in range (geo.Nx):
        for j in range(1, geo.Ny-1):

            
            Fields.v_psu[i,j] = 0.7*Fields.v_psu[i,j] + (1 - 0.7)*Fields.v_old[i,j]

            Fields.v[i,j] = Fields.v_psu[i,j] + (Coupler.d_ns[i,j-1]*(Fields.P_prime[i,j]-Fields.P_prime[i,j+1])) 



    print('v final results')
    print(np.array2string(np.flipud(Fields.v.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('u final results')
    print(np.array2string(np.flipud(Fields.u.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))
    print('P final results')
    print(np.array2string(np.flipud(Fields.P.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))

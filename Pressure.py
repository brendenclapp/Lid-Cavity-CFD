import numpy as np
from scipy.linalg import solve_banded

def Coeff_P(geo, var, Faces, Fields, Coupler):

    Nx = geo.Nx
    Ny = geo.Ny

    ip = 0
    jp = 0


    for ip in range (Nx):
        for jp in range (Ny):

            if Faces.a_P_u[ip,jp] < 1e-12:

                Coupler.d_we[ip,jp] = 0

            else:

                Coupler.d_we[ip,jp] = geo.dy / Faces.a_P_u[ip,jp]


    for ipp in range (Nx):
        for jpp in range (Ny):

            if Faces.a_P_v[ipp,jpp] < 1e-12:

                 Coupler.d_ns[ipp,jpp] = 0

            else:

                Coupler.d_ns[ipp,jpp] = geo.dx / Faces.a_P_v[ipp,jpp]
    
    Coupler.a_NS_P[:,0] = 0
    Coupler.a_NS_P[:,-1] = 0

    i = 0
    j = 0
    for i in range (Nx+1):
        for j in range (Ny):

            Coupler.a_EW_P[i,j] = var.rho * Coupler.d_we[i,j] * geo.dy

    i = 0
    j = 0
    for i in range (Nx):
        for j in range (Ny+1):

            Coupler.a_NS_P[i,j] = var.rho * Coupler.d_ns[i,j] * geo.dx

    print('EW')
    print(np.array2string(np.flipud(Coupler.a_EW_P.T),formatter={'float_kind': lambda x: f"{x:10.3f}"}, max_line_width= 1000000000))
    print('NS')
    print(np.array2string(np.flipud(Coupler.a_NS_P.T),formatter={'float_kind': lambda x: f"{x:10.3f}"}, max_line_width= 1000000000))

    i = 0
    j = 0
    for i in range(Nx):
        for j in range(Ny):

            Coupler.b[i,j] = var.rho*(((Fields.u_psu[i,j]*geo.dy) - (Fields.u_psu[i+1,j]*geo.dy)) + ((Fields.v_psu[i,j]*geo.dx) - (Fields.v_psu[i,j+1]*geo.dx)))
            
    print('b')
    print(np.array2string(np.flipud(Coupler.b.T),formatter={'float_kind': lambda x: f"{x:10.3f}"}, max_line_width= 1000000000))

    i = 0
    j = 0
    for i in range (Nx):
        for j in range (Ny):

            Coupler.a_P_P[i,j] = Coupler.a_EW_P[i,j] + Coupler.a_EW_P[i+1,j] + Coupler.a_NS_P[i,j] + Coupler.a_NS_P[i,j+1]


def TDMA_P(geo, var, Faces, Fields, Coupler, Tri):

    ip = 0
    jp = 0

    Fields.P_prime_old = np.copy(Fields.P_prime) 
    for ip in range (geo.Nx):
        for jp in range (geo.Ny):


            Tri.upper_P[jp] = -Coupler.a_NS_P[ip,jp+1]
            Tri.diag_P[jp] = Coupler.a_P_P[ip,jp]
            Tri.lower_P[jp] = -Coupler.a_NS_P[ip,jp]
            Tri.RHS_P[jp] = Coupler.b[ip,jp]

            #West Zero-Gradiant
            if ip == 0:

                Tri.diag_P[jp] -= Coupler.a_EW_P[ip,jp]
                Tri.RHS_P[jp] += Coupler.a_EW_P[ip+1,jp] * Fields.P_prime_old[ip+1,jp]

            #East Dirchlet anchor
            elif ip == geo.Nx-1:

                Tri.RHS_P[jp] += Coupler.a_EW_P[ip+1,jp] * 0
                Tri.RHS_P[jp] += Coupler.a_EW_P[ip,jp] * Fields.P_prime[ip-1,jp]

            #Interior
            else:

                Tri.RHS_P[jp] += Coupler.a_EW_P[ip,jp]* Fields.P_prime[ip-1,jp]
                Tri.RHS_P[jp] += Coupler.a_EW_P[ip+1,jp]*Fields.P_prime_old[ip+1,jp] 

            #South wall Zero-Gradiant
            if jp == 0:

                Tri.lower_P[jp] = 0
                Tri.diag_P[jp] -= Coupler.a_NS_P[ip,jp]

            #North wall Zero-Gradiant
            elif jp == geo.Ny-1:

                Tri.upper_P[jp] = 0
                Tri.diag_P[jp] -= Coupler.a_NS_P[ip,jp+1]

        ab = np.zeros((3, geo.Ny))
        ab[0, 1:] = Tri.upper_P[:-1]
        ab[1,:] = Tri.diag_P
        ab[2,:-1] = Tri.lower_P[1:]
        sol = solve_banded((1,1), ab, Tri.RHS_P)
        Fields.P_prime[ip,:] = sol
            

    print('P_TDMA results')
    print(np.array2string(np.flipud(Fields.P_prime.T),formatter={'float_kind': lambda x: f"{x:8.3f}"}, max_line_width= 1000000000))



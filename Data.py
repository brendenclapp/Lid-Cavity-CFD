import numpy as np


class Geometry:

    def __init__ (self):

        self.Nx = 80             # number of cells in x direction
        self.Ny = 80             # number of cells in y direction

        self.lx = 0.01           # (m) length of computational domain
        self.ly = 0.01           # (m) height of computational domain

        self.dx = self.lx / ( self.Nx - 1 )     # (m) Length of cell
        self.dy = self.ly / ( self.Ny - 1 )     # (m) Height of cell

        self.ar = self.dx / self.dy             # Calculation Simplification
        self.ra = 1/self.ar                     # Calculation Simplification


class Variables:

    def __init__(self, geo):

        
        self.rho = 1000                        # (kg/m^3) density
        self.mu = 0.001                        # (Pa*s) dynamic viscosity

        self.u_lid = 0.01                      # (m/s) speed of cavity lid                  

        self.res_u = 1                         # residual of u solution (1 is placeholder)
        self.res_v = 1                         # residual of v solution (1 is placeholder)
        self.res_P = 1                         # residual of P solution (1 is placeholder)

        self.utol = 1e-6                       # u solution tolerance
        self.vtol = 1e-6                       # v solution tolerance
        self.ptol = 6.5e-5                     # P solution tolerance

        self.itt = 1                           # Outer loop itteration counter

class Fields:        

    def __init__(self, geo):

        self.P = np.full((geo.Nx, geo.Ny), 0.0)             # Pressure field
        self.u = np.full((geo.Nx+1, geo.Ny), 0.0)           # X-momentum field
        self.v = np.full((geo.Nx, geo.Ny+1), 0.0)           # Y-momentum field
        self.u_tilde = np.zeros((geo.Nx+1, geo.Ny))         # X-momentum correction field
        self.v_tilde = np.zeros((geo.Nx, geo.Ny+1))         # Y-momentum correction field
        self.Pp = np.zeros((geo.Nx, geo.Ny))                # Pressure correction field


class Faces:

    def __init__(self, geo):


        self.a_w_u = np.zeros((geo.Nx+1, geo.Ny))           # western variable coefficient of u CV for U solver
        self.a_e_u = np.zeros_like(self.a_w_u)              # eastern variable coefficient of u CV for U solver
        self.a_n_u = np.zeros_like(self.a_w_u)              # northern variable coefficient of u CV for U solver
        self.a_s_u = np.zeros_like(self.a_w_u)              # southern variable coefficient of u CV for U solver
        self.a_P_u = np.zeros_like(self.a_w_u)              # center variable coefficient of u CV for U solver

        self.dPdx = np.zeros_like(self.a_w_u)               # horizontal pressure gradiant 

        self.a_w_v = np.zeros((geo.Nx, geo.Ny+1))           # western variable coefficient of v CV for V solver
        self.a_e_v = np.zeros_like((self.a_w_v))            # eastern variable coefficient of v CV for V solver
        self.a_n_v = np.zeros_like(self.a_w_v)              # northern variable coefficient of v CV for V solver
        self.a_s_v = np.zeros_like(self.a_w_v)              # southern variable coefficient of v CV for V solver
        self.a_P_v = np.zeros_like(self.a_w_v)              # center variable coefficient of v CV for V solver

        self.dPdy = np.zeros_like(self.a_w_v)               # vertical pressure gradiant 

        self.a_w_P = np.zeros((geo.Nx, geo.Ny))             # western variable coefficient of P CV for P solver
        self.a_e_P = np.zeros((geo.Nx, geo.Ny))             # eastern variable coefficient of P CV for P solver
        self.a_n_P = np.zeros((geo.Nx, geo.Ny))             # northern variable coefficient of P CV for P solver
        self.a_s_P = np.zeros((geo.Nx, geo.Ny))             # southern variable coefficient of P CV for P solver
        self.a_P_P = np.zeros((geo.Nx, geo.Ny))             # center variable coefficient of P CV for P solver

        self.b = np.zeros((geo.Nx, geo.Ny))

class TDMA:

    def __init__ (self,geo):

        # TDMA column arrays for the u,v and P ADI algorithms

        # X-momentum
        self.upper_col = np.zeros((geo.Ny))
        self.diag_col = np.zeros((geo.Ny))
        self.lower_col = np.zeros((geo.Ny))
        self.RHS_col = np.zeros((geo.Ny))

        self.upper_row = np.zeros((geo.Ny-1))
        self.diag_row = np.zeros((geo.Ny-1))
        self.lower_row = np.zeros((geo.Ny-1))
        self.RHS_row = np.zeros((geo.Ny-1))

        # Y-momentum
        self.upper_colv = np.zeros((geo.Ny-1))
        self.diag_colv = np.zeros((geo.Ny-1))
        self.lower_colv = np.zeros((geo.Ny-1))
        self.RHS_colv = np.zeros((geo.Ny-1))

        self.upper_rowv = np.zeros((geo.Ny))
        self.diag_rowv = np.zeros((geo.Ny))
        self.lower_rowv = np.zeros((geo.Ny))
        self.RHS_rowv = np.zeros((geo.Ny))

        # Pressure
        self.upper_P = np.zeros((geo.Ny))
        self.diag_P = np.zeros((geo.Ny))
        self.lower_P = np.zeros((geo.Ny))
        self.RHS_P = np.zeros((geo.Ny))


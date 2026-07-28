import numpy as np


class Geometry:

    def __init__ (self):

        self.Nx = 8
        self.Ny = 8

        self.lx = 0.01            # meters
        self.ly = 0.01           # meters

        self.dx = self.lx / ( self.Nx - 1 )
        self.dy = self.ly / ( self.Ny - 1 )

        self.ar = self.dx / self.dy
        self.ra = 1/self.ar


class Variables:

    def __init__(self, geo):

        
        self.rho = 1000
        self.mu = 0.001

        self.u_lid = 0.01     # m/s
        self.v_inlet = 0.0

        self.De = (self.mu / geo.dx) * geo.dy
        self.Dw = (self.mu / geo.dx) * geo.dy
        self.Dn = (self.mu / geo.dy) * geo.dx
        self.Ds = (self.mu / geo.dy) * geo.dx

class Fields:        

    def __init__(self, geo):

        self.P = np.full((geo.Nx, geo.Ny), 0.0)
        self.u = np.full((geo.Nx+1, geo.Ny), 0.0)
        self.v = np.full((geo.Nx, geo.Ny+1), 0.0)
        self.u_psu = np.zeros((geo.Nx+1, geo.Ny))              # field that holds psuedo u, the TDMA u solution, which has not yet been pressure corrected
        self.u_tilde = np.zeros((geo.Nx+1, geo.Ny))
        self.v_tilde = np.zeros((geo.Nx, geo.Ny+1))
        self.u_old = np.zeros((geo.Nx+1, geo.Ny))
        self.v_old = np.zeros((geo.Nx, geo.Ny+1))
        self.v_psu = np.zeros((geo.Nx, geo.Ny+1))
        self.Pp = np.zeros((geo.Nx, geo.Ny))
        self.P_prime_old = np.zeros((geo.Nx, geo.Ny))


class Faces:

    def __init__(self, geo):


        self.F_we = np.zeros((geo.Nx+1, geo.Ny))           # West East faces for the scalar CV
        self.F_ns = np.zeros((geo.Nx, geo.Ny+1))           # North South faces for the scalar CV
        self.F_we_psu = np.zeros((geo.Nx+1, geo.Ny))           # West East faces for the scalar CV
        self.F_ns_psu = np.zeros((geo.Nx, geo.Ny+1))           # North South faces for the scalar CV

        self.Fw_u = np.zeros((geo.Nx+1, geo.Ny))            # Convection Flux for western u CV 
        self.Fe_u = np.zeros((geo.Nx+1, geo.Ny))            # Convection Flux for eastern u CV 
        self.Fn_u = np.zeros((geo.Nx+1, geo.Ny))            # Convection Flux for northern u CV 
        self.Fs_u = np.zeros((geo.Nx+1, geo.Ny))            # Convection Flux for southern u CV 

        self.a_w_u = np.zeros((geo.Nx+1, geo.Ny))           # western variable coefficient of u CV for U solver
        self.a_e_u = np.zeros_like(self.a_w_u)              # eastern variable coefficient of u CV for U solver
        self.a_n_u = np.zeros_like(self.a_w_u)              # northern variable coefficient of u CV for U solver
        self.a_s_u = np.zeros_like(self.a_w_u)              # southern variable coefficient of u CV for U solver
        self.a_P_u = np.zeros_like(self.a_w_u)              # center variable coefficient of u CV for U solver

        self.dPdx = np.zeros_like(self.a_w_u)               # horizontal pressure gradiant located at each ( shared w/ ) U node

        self.Fw_v = np.zeros((geo.Nx, geo.Ny+1))            # Convection Flux for western v CV 
        self.Fe_v = np.zeros((geo.Nx, geo.Ny+1))            # Convection Flux for eastern v CV 
        self.Fn_v = np.zeros((geo.Nx, geo.Ny+1))            # Convection Flux for northern v CV 
        self.Fs_v = np.zeros((geo.Nx, geo.Ny+1))            # Convection Flux for southern v CV 

        self.a_w_v = np.zeros((geo.Nx, geo.Ny+1))           # western variable coefficient of v CV for V solver
        self.a_e_v = np.zeros_like((self.a_w_v))            # eastern variable coefficient of v CV for V solver
        self.a_n_v = np.zeros_like(self.a_w_v)              # northern variable coefficient of v CV for V solver
        self.a_s_v = np.zeros_like(self.a_w_v)              # southern variable coefficient of v CV for V solver
        self.a_P_v = np.zeros_like(self.a_w_v)              # center variable coefficient of v CV for V solver

        self.dPdy = np.zeros_like(self.a_w_v)               # vertical pressure gradiant located at each ( shared w/ ) V node

        self.a_w_P = np.zeros((geo.Nx, geo.Ny))          
        self.a_e_P = np.zeros((geo.Nx, geo.Ny))          
        self.a_n_P = np.zeros((geo.Nx, geo.Ny))            
        self.a_s_P = np.zeros((geo.Nx, geo.Ny))           
        self.a_P_P = np.zeros((geo.Nx, geo.Ny))            

class TDMA:

    def __init__ (self,geo):

        self.upper_col = np.zeros((geo.Ny))
        self.diag_col = np.zeros((geo.Ny))
        self.lower_col = np.zeros((geo.Ny))
        self.RHS_col = np.zeros((geo.Ny))

        self.upper_row = np.zeros((geo.Ny-1))
        self.diag_row = np.zeros((geo.Ny-1))
        self.lower_row = np.zeros((geo.Ny-1))
        self.RHS_row = np.zeros((geo.Ny-1))

        self.upper_colv = np.zeros((geo.Ny-1))
        self.diag_colv = np.zeros((geo.Ny-1))
        self.lower_colv = np.zeros((geo.Ny-1))
        self.RHS_colv = np.zeros((geo.Ny-1))

        self.upper_rowv = np.zeros((geo.Ny))
        self.diag_rowv = np.zeros((geo.Ny))
        self.lower_rowv = np.zeros((geo.Ny))
        self.RHS_rowv = np.zeros((geo.Ny))

        self.upper_v = np.zeros((geo.Ny-1))
        self.diag_v = np.zeros((geo.Ny-1))
        self.lower_v = np.zeros((geo.Ny-1))
        self.RHS_v = np.zeros((geo.Ny-1))

        self.upper_P = np.zeros((geo.Ny))
        self.diag_P = np.zeros((geo.Ny))
        self.lower_P = np.zeros((geo.Ny))
        self.RHS_P = np.zeros((geo.Ny))


class Coupler:

    def __init__(self,geo):

        self.d_we = np.zeros((geo.Nx+1, geo.Ny))
        self.d_ns = np.zeros((geo.Nx, geo.Ny+1))

        self.a_EW_P = np.zeros((geo.Nx+1, geo.Ny))
        self.a_NS_P = np.zeros((geo.Nx, geo.Ny+1))
        self.a_P_P = np.zeros((geo.Nx, geo.Ny))

        self.b = np.zeros((geo.Nx, geo.Ny))
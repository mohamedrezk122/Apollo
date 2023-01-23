from reader import *
import numpy as np 
from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
import json 

# Frame of Reference ECI 
# units in km , kg, and s
cos = np.cos
sin = np.sin 
pi  = np.pi 

def rad2deg(angle):
    return (180*angle)/pi

def deg2rad(angle):
    return (pi*angle)/180

def cosd(angle):
    return np.cos(deg2rad(angle))
def sind(angle):
    return np.sin(deg2rad(angle))

# print(cosd(30))
class Propagator:

    def __init__(self , r , v ):

        self.mu = 6.6743e-11 * 5.972e24 * 1e-9  # GM in km^3 s^-2
        self.r = r  # position vector
        self.v = v  # velocity vector
        self.r_mag = linalg.norm(self.r) # magnitude of the position vector
        self.v_mag = linalg.norm(self.v) # magnitude of the velocity vector
        self.x     = np.array([1,0,0])   # x-axis (vernal equinox) 
        self.y     = np.array([0,1,0])   # cross product of x and z
        self.z     = np.array([0,0,1])   # z-axis (north pole) 
        self.h     = np.cross(self.r , self.v) # angular momentum vector
        self.n     = np.cross(self.z,self.h)  # line of nodes vector
        self.e     = self.Compute_eccentricity()    # eccentricity vector
        self.e_mag = linalg.norm(self.e)   # magnitude of eccentricity vector
        self.h_mag = linalg.norm(self.h)   # magnitude of angular momentum vector
        self.n_mag = linalg.norm(self.n)
        self.i     = self.Compute_inclination()# orbit inclination
        self.E     = (0.5*self.v_mag**2) - (self.mu/self.r_mag) # motion energy
        self.a     = - self.mu/(2*self.E)
        self.mean_motion = np.sqrt(self.mu/self.a**3)
        self.omega = self.Compute_RAAN() #right ascension of ascending node
        self.perigee= self.Compute_arg_perigee()
        self.f     = self.Compute_true_anomaly() 
        self.period = 2*pi *np.sqrt(self.a**3/self.mu) # in seconds


    def Compute_eccentricity(self):
        return (1/self.mu)* (np.multiply(self.v_mag**2 -(self.mu/self.r_mag) , self.r) - \
                            np.multiply(np.dot(self.r,self.v),self.v))

    def Compute_inclination(self):
        return rad2deg(np.arccos(np.dot(np.multiply((1/self.h_mag),self.h),self.z))) 

    def Compute_RAAN(self):
        omega = rad2deg(np.arccos(np.dot(np.multiply((1/self.n_mag),self.n),self.x))) 
        # quadrant check 
        if self.n[1] < 0 :
            omega = 360 - omega    
        return omega

    def Compute_arg_perigee(self):
        omg = rad2deg(np.arccos(np.dot(self.n,self.e)/(self.n_mag*self.e_mag)))
        # quadrant check 

        if np.dot(self.e, self.z) < 0 :
            omg = 360 - omg    
        return omg
    def Compute_true_anomaly(self):
        f = rad2deg(np.arccos(np.dot(self.r,self.e)/(self.r_mag*self.e_mag)))
        # quadrant check 
        if np.dot(self.r,self.v) < 0:
            f = 360 - f    
        return f

    def Compute_keplers(self, delta_t):
        # print((1-self.e_mag)/(1+self.e_mag))
        E_t1 = 2*np.arctan((np.sqrt((1-self.e_mag)/(1+self.e_mag))*np.tan(deg2rad(self.f)/2)))
        M_t1 = E_t1 - self.e_mag*sin(E_t1)
        M_t2 = M_t1 + self.mean_motion*(delta_t)

        # Netwon-Raphson iterative approach 
        max_iter = 50 
        i = 0
        E_k = M_t2 
        while abs(M_t2-E_k+self.e_mag*sin(E_k)) > 0.001 and i < max_iter :
            E_k = E_k - (M_t2-E_k +self.e_mag*sin(E_k))/(self.e_mag*cos(E_k) - 1)
            i+=1
        E_t2 = E_k
        f_t2 = rad2deg(2*np.arctan((np.sqrt((1+self.e_mag)/(1-self.e_mag))*np.tan(E_t2/2))))
        r_new = np.array([
            self.r_mag*(cosd(self.omega)*cosd(self.perigee+f_t2)-sind(self.omega)*sind(self.perigee+f_t2)*cosd(self.i)),
            self.r_mag*(sind(self.omega)*cosd(self.perigee+f_t2)+cosd(self.omega)*sind(self.perigee+f_t2)*cosd(self.i)),
            self.r_mag*(sind(self.perigee+f_t2)*sind(self.i)) ])

        
        v_new = np.array([ 
                (-self.mu/self.h_mag)*(cosd(self.omega)*(sind(self.perigee+f_t2)+\
                self.e_mag*sind(self.perigee))+ \
                sind(self.omega)*(cosd(self.perigee+f_t2)+\
                self.e_mag*cosd(self.perigee))*cosd(self.i)),

            (-self.mu/self.h_mag)*(sind(self.omega)*(sind(self.perigee+f_t2)+\
                self.e_mag*sind(self.perigee))- \
                cosd(self.omega)*(cosd(self.perigee+f_t2)+\
                self.e_mag*cosd(self.perigee))*cosd(self.i)),

            (self.mu/self.h_mag)*(cosd(self.perigee+f_t2)+\
                self.e_mag*cosd(self.perigee)*sind(self.i))
            ])
        return r_new ,v_new

    def Rx(self,a):
        # rotation matrix along x-axis
        return np.array([[1     ,   0    ,    0   ],
                        [0     , cos(a) , -sin(a) ],
                        [0     , sin(a) ,  cos(a) ]])

    def Ry(self,a):
        # rotation matrix along y-axis
        return np.array([[cos(a), 0     ,  sin(a) ],
                        [0     , 1      ,    0    ],
                        [-sin(a), 0      , cos(a) ]])

    def Rz(self,a):
        # rotation matrix along z-axis
        return np.array([[cos(a), -sin(a),    0   ],
                        [sin(a), cos(a) ,    0    ],
                        [0     ,    0   ,    1    ]])


# r ,v = get_state_vector(2022,9,29,20,14,2)

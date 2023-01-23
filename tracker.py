import numpy as np 
import matplotlib.pyplot as plt
import math as m
from matplotlib.animation import FuncAnimation 
import pyproj
cos = np.cos
sin = np.sin 
pi  = np.pi 

def rad2deg(angle):
    return (180*angle)/pi

def deg2rad(angle):
    return (pi*angle)/180

class Tracker:

    def __init__(self,y,m,d, ISS_position , ISS_ascension, ISS_inclination,ISS_perigee):

        self.P = ISS_position
        self.omega = ISS_perigee
        self.Omega = ISS_ascension
        self.i = ISS_inclination
        self.jd  = self.JD(y,m,d)
        self.tau = (self.jd-2451545)/365250 # julian millennia from epoch j2000
        self.T   =  self.tau *10 # julian centuries from epoch j2000
        lonTerms,latTerms,radTerms = self.readVSOP()
        hL, hB, self.R = self.HeliocentricLBR(lonTerms,latTerms,radTerms)
        self.rad = 1.496e+8*self.R
        self.theta , self.beta = self.GeocentricLB(hL,hB) 
        self.FK5Correction()
        self.NutationCorrection()
        self.lamb = self.AberrationCorrection()
        self.alpha , self.delta = self.EclipticToEquatorial()
        self.Li = self.Compute_Li()
        self.Lo = self.Compute_Lo()
        
        # print("radius = " , self.rad)
        # print("alpha =" , self.alpha)
        # print('delta = ' , self.delta)
        # print("Li =" , self.Li)
        # print("Lo =" , self.Lo)

    def JD(self,y,m,d):
        # julian days from epoch j2000
        # m in numbers
        if m in [1,2]:
            y = y - 1
            m = m + 12 
        A = int(y/100)
        B = 2 - A + int(A/4)
        jd = int(365.25*(y+4716)) + int(30.6001*(m+1)) \
                + d + B - 1524.5
        return jd

    def readVSOP(self):

        # Read the periodic terms for a heliocentric ecliptical planet position from a VSOP87D.* file.
        
        fileName = './VSOP87D.ear'
        inFile = open(fileName,'r')
        
        import fortranformat as ff
        formatHeader = ff.FortranRecordReader('(40x,I3, 16x,I1,I8)')         # Block header format
        formatBody   = ff.FortranRecordReader('(79x,F18.11,F14.11,F20.11)')  # Block body format
        
        lonTerms=[]; latTerms=[]; radTerms=[]
        
        for iBlock in range(3*6):  # 3 variables (l,b,r), up to 6 powers (0-5)
            line = inFile.readline()
            var,power,nTerm = formatHeader.read(line)
            if line == '': break  # EoF
            
            for iLine in range(nTerm):
                line = inFile.readline()
                a,b,c = formatBody.read(line)
                if var == 1: lonTerms.append([power, a,b,c])  # var=1: ecliptic longitude
                if var == 2: latTerms.append([power, a,b,c])  # var=2: ecliptic latitude
                if var == 3: radTerms.append([power, a,b,c])  # var=3: radial distance
                
        return lonTerms,latTerms,radTerms

    def HeliocentricLBR(self,lonTerms,latTerms,radTerms):

        tau = self.tau
        L = 0.0 ; B = 0.0 ;  R = 0.0 

        for term in lonTerms:
            costerm = term[1]*cos(term[2]+term[3]*tau)
            L += costerm * tau**term[0]

        for term in latTerms:
            costerm = term[1]*cos(term[2]+term[3]*tau)
            B += costerm * tau**term[0]

        for term in radTerms:
            costerm = term[1]*cos(term[2]+term[3]*tau)
            R += costerm * tau**term[0]

        return rad2deg(L%(2*pi)),rad2deg(B),R

    def GeocentricLB(self,L,B):

        # convert heliocentric location of (sun to earth) to 
        # geocentric location of (earth to sun)
        return L+180 , -B

    def FK5Correction(self):
        gamma = deg2rad(self.theta - 1.397*self.T-0.00031*self.T**2)
        delta_theta = -0.09033/3600
        delta_beta = (0.03916/3600) * (cos(gamma)-sin(gamma))
        self.theta = self.theta + delta_theta
        self.beta  = self.beta  + delta_beta 

    def NutationCorrection(self):
        m = deg2rad(280.4665 + 36000.7698*self.T)
        mm= deg2rad(218.3165 + 481267.8813*self.T)
        omegar = 125.04452 - 1934.136261*self.T + 0.0020708*self.T**2\
                    +(self.T**3)/450000
        delta_epsi = (-17.2*sin(omegar)-1.32*sin(2*m)\
                    -0.23*sin(2*mm)+0.21*sin(2*omegar))/3600
        self.theta = self.theta + delta_epsi
        return m ,mm, omegar

    def AberrationCorrection(self):
        lamb = deg2rad(self.theta - (20.4898/3600)/(self.R))
        return lamb 

    def EclipticToEquatorial(self):
        m , mm , omegar = self.NutationCorrection()
        u = self.T/100
        delta_epsilon =  (9.2*cos(omegar)+0.57*cos(2*m)\
                    +0.1*cos(2*mm)-0.09*cos(2*omegar))/3600
        epsilon_0 = 23.4392911 - (4680.93*u -1.55*u**2 -\
                    1999.25*u**3 -51.38*u**4 -249.67*u**5 -39.05*u**6 + \
                    7.12*u**7 + 27.87*u**8 + 5.79*u**9 + 2.45*u**10)/3600 
        epsilon = deg2rad(epsilon_0 + delta_epsilon)

        ascension = np.arctan((sin(self.lamb)*cos(epsilon)\
                 - np.tan(deg2rad(self.beta))* sin(epsilon) )/cos(self.lamb))

        declination = np.arcsin((sin(deg2rad(self.beta))*cos(epsilon))+\
                            (cos(deg2rad(self.beta))*sin(epsilon)*sin(self.lamb)))
        return ascension , declination

    def Compute_Li(self):
        # a vector points from earth center to sun  
        mat = np.array([cos(self.alpha)*cos(self.delta) ,
                        cos(self.delta)*sin(self.alpha) ,
                        sin(self.delta) ])
        L_i = np.multiply(self.rad , mat)

        return L_i
    def Compute_Lo(self):
         # a vector points from satellite center to sun 

        a = deg2rad(self.i-90) 
        b = deg2rad(-90-self.omega)
        c = deg2rad(self.Omega)

        R_x = np.array([[1     ,   0    ,    0    ],
                        [0     , cos(a) , -sin(a) ],
                        [0     , sin(a) ,  cos(a) ]])

        R_y = np.array([[cos(b), 0      ,  sin(b) ],
                        [0     , 1      ,    0    ],
                        [-sin(b), 0      , cos(b) ]])

        R_z = np.array([[cos(c), -sin(c),    0    ],
                        [sin(c), cos(c) ,    0    ],
                        [0     ,    0   ,    1    ]])

        R = np.matmul(R_y,R_x,R_z)
        Lo = np.matmul(R,self.Li)+self.P
        return Lo



# import json
# import urllib.request as urllib2
# req = urllib2.Request("http://api.open-notify.org/iss-now.json")
# response = urllib2.urlopen(req)

# obj = json.loads(response.read())
# print (obj['iss_position']['latitude'], obj['iss_position']['longitude'])

# url ='http://celestrak.org/NORAD/elements/stations.txt'
# req = urllib2.Request(url)
# response = urllib2.urlopen(req)
# data = response.read()
# tle = data.decode().split("\n")[0:3]
# line2 = tle[2]

# perigee     = float(line2[34:42])
# inclination = float(line2[8:16])
# ISS_ascension = float(line2[17:25])
# ecc = line2[26:33]


# p = np.array([5673.170960667850 , 3295.164548032450 , -1775.016192577940 ])
# ts= TrackSun(2022,9,2 , p,ISS_ascension,inclination,perigee )

# print(p)
# """
# - write an xml file handler to get the live state vector
# - plot and test
# """

# figure = plt.figure()
#     # figure.set_tight_layout(True)   
# ax =figure.add_subplot(111 , projection='3d')
# #earth

# # print(np.sqrt(ts.Li[0]**2 + ts.Li[1]**2 + ts.Li[2]**2))
# ax.scatter3D(0,0,0,color='blue')
# ax.scatter3D(p[0],p[1],p[2],color ="green")
# ax.scatter3D(ts.Li[0],ts.Li[1],ts.Li[2],color ="orange",s=100)
# ax.plot([0,p[0]],[0,p[1]],[0,p[2]],color ='yellow')
# ax.plot([ts.Li[0],ts.Lo[2]],[ts.Li[1],ts.Lo[2]],[ts.Li[2],ts.Lo[2]] , color='red')
# ax.plot([ts.Li[0],0],[ts.Li[1],0],[ts.Li[2],0] , color='green')
# plt.show()



# X, Y, Z = np.meshgrid(np.arange(-0.8, 1, 0.2),
#                       np.arange(-0.8, 1, 0.2),
#                       np.arange(-0.8, 1, 0.8))

# U = np.sin(np.pi * X) * np.cos(np.pi * Y) * np.cos(np.pi * Z)
# V = -np.cos(np.pi * X) * np.sin(np.pi * Y) * np.cos(np.pi * Z)
# W = (np.sqrt(2.0 / 3.0) * np.cos(np.pi * X) * np.cos(np.pi * Y) *
#      np.sin(np.pi * Z))

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.set_title("3D Quiver plot")

# #drawing quiver plot
# ax.quiver(X, Y, Z, U, V, W, length=0.1, normalize=True)

# plt.show()


# x =  ((MAP_WIDTH/360.0) * (180 + lon))
# y = ((MAP_HEIGHT/180.0) * (90 - lat))

# pastx = []
# pasty = []
# pastz = []
    



  

# def get():
#     req = urllib2.Request("http://api.open-notify.org/iss-now.json")
#     response = urllib2.urlopen(req)
#     obj = json.loads(response.read())
#     lat ,  lon = deg2rad(float(obj['iss_position']['latitude'])), deg2rad(float(obj['iss_position']['longitude']))
#     R = 6375
#     x = R * np.cos(lat) * np.cos(lon)
#     y = R * np.cos(lat) * np.sin(lon)
#     z = R *sin(lat)
#     print("me: ", x,y,z)
#     # hae = 408  # altitude
#     # transformer = pyproj.Transformer.from_crs(
#     # {"proj":'latlong', "ellps":'WGS84', "datum":'WGS84'},
#     # {"proj":'geocent', "ellps":'WGS84', "datum":'WGS84'},
#     # )
#     # x ,y, z = transformer.transform(lon,lat,hae,radians = False)
#     # print("lib: ", x/1000,y/1000,z/1000)
#     # print("-"*50)
#     return x,y,z


# def frames():
#     while True:
#         yield get()


# x = []
# y = []
# z = []
# def animate(args):

#     pastx.append(args[0])
#     pasty.append(args[1])
#     pastz.append(args[2])
    
#     return ax.plot(pastx,pasty,pastz,color='blue')


# anim = FuncAnimation(figure, animate, frames=frames, interval=1000)
# plt.show()
from propagator import Propagator
from tracker import Tracker
from reader import * 
import numpy as np 
from numpy import linalg
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation 
import json 

# y,m,d, ISS_position , ISS_ascension, ISS_inclination,ISS_perigee
figure = plt.figure()
ax =figure.add_subplot(111)
R = 6375
us, vs = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
xs = np.cos(us)*np.sin(vs)*R
ys = np.sin(us)*np.sin(vs)*R
zs = np.cos(vs)*R

ax.set_aspect('auto')

pastx = []
pasty = []
pastz = []
t = 0 
dt  =  50
y = 2022
m = 9
d  =  30
h = 0
mi = 2
_,r,v= get_state_vector(y,m,d,h,mi,3)

i = 0 

def align_vectors(a, b):
    b = b / np.linalg.norm(b) # normalize a
    a = a / np.linalg.norm(a) # normalize b
    v = np.cross(a, b)
    c = np.dot(a, b)

    v1, v2, v3 = v
    h = 1 / (1 + c)

    Vmat = np.array([[0, -v3, v2],
                  [v3, 0, -v1],
                  [-v2, v1, 0]])

    R = np.eye(3, dtype=np.float64) + Vmat + (Vmat.dot(Vmat) * h)
    return R
    
def frames():
    global t
    while True:
        t += dt
        g = t /(3600*24)
        P = Propagator(r,v)
        r1,v1 = P.Compute_keplers(t)
        track = Tracker(y,m,d+g,r1,P.omega,P.i , P.perigee)

        yield r1 , track.Li, track.Lo 

def animate(args):


    r2 = args[0]
    Li = args[1]
    Lo = args[2]
    # print("Lo" ,Lo)
    pastx.append(r2[0])
    pasty.append(r2[1])
    pastz.append(r2[2])
    ax.cla() 
    # ax.set_xlim(-7000,7000)
    # ax.set_xlim(-7000,7000)

    # ax.plot_wireframe(xs, ys, zs, color="r")
    print([r2[0], Li[0]])
    ax.plot(pastx,pasty,color='blue')
    ax.plot([0, Li[0]], [0, Li[1]] , color ="green")
    # ax.plot([, Li[0]], [0, Li[1]],[0, Li[2]] , color ="green")
    ax.scatter(0,0 , s=80 , color="blue")
    ax.scatter(Li[0],Li[1] , s=500 , color='yellow')

    # ax.scatter3D(L)
    return ax.plot([r2[0], Li[0]], [r2[1], Li[1]] , color="orange")





anim = FuncAnimation(figure, animate, frames=frames, interval=1000)
plt.show()
import numpy as np 
import matplotlib.pyplot as plt


figure = plt.figure()
ax =figure.add_subplot(111, projection="3d")

r1 = np.array([2,3,4])
r2 = np.array([1,6,7])

cos  = np.cos
sin =  np.sin 

ax.plot([0,r1[0]],[0,r1[1]] , [0,r1[2]] , color="green")

ax.plot([0,r2[0]],[0,r2[1]] , [0,r2[2]] , color="blue")


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


R = align_vectors(r1,r2)

r3 = R@r1

ax.plot([0,r3[0]],[0,r3[1]] , [0,r3[2]] , color="red")
plt.show()


import numpy as np
import os
import matplotlib.pyplot as plt

nStep = 300
iMax = 100
dt = 0.001
dx = 0.01

iMax = iMax + 2 #including ghost cell

gamma = 1.403

u = np.zeros(iMax)
rho = np.zeros(iMax)
p = np.zeros(iMax)
e = np.zeros(iMax)
x = np.zeros(iMax)

lBdry = np.zeros(3) # boundary ghost cells
rBdry = np.zeros(3) # boundary ghost cells
qf = np.zeros((iMax,3))
Qc = np.zeros((iMax,3))

def initialCondition():
    for i in range(iMax):
        if i <= iMax*0.5:
            rho[i] = 1.0
            p[i] = 1.0
        else:
            rho[i] = 0.125
            p[i] = 0.1
        e[i] = p[i]/(gamma-1)+rho[i]*(u[i]**2)/2
        x[i] = i*dx-dx/2
    
    lBdry[0] = rho[0]
    rBdry[0] = rho[iMax-1]
    lBdry[1] = u[0]*rho[0]
    rBdry[1] = u[iMax-1]*rho[iMax-1]
    lBdry[2] = e[2]
    rBdry[2] = e[iMax-1]

    for i in range(iMax):
        qf[i][0] = u[i]
        Qc[i][0] = rho[i]
        qf[i][1] = rho[i]
        Qc[i][1] = u[i]*rho[i]
        qf[i][2] = p[i]
        Qc[i][2] = e[i]

def boundaryCondition():
    global Qc
    for i in range(3):
        Qc[0][i] = 2*lBdry[i]-Qc[1][i] # left boundary
        Qc[iMax-1][i] = Qc[iMax-2][i] # right boundary

def calQ():
    global Qc
    calRes()

    Qc = np.array(Qc)
    lo_R = np.array(Res)

    for i in range(1,iMax-1):
        Qc[i] = Qc[i] -dt/dx*lo_R[i]
    
    Qc.tolist()
    boundaryCondition()

def calRes():
    global Res
    Res = np.zeros((iMax,3))
    fvs()
    for i in range(1,iMax-1):
        Res[i] = fluxPlus[i] - fluxPlus[i-1]
    Res.tolist()

def fvs():
    global fluxPlus
    fluxPlus = np.zeros((iMax,3)) # fluxPlus[i] is numerical flux between cell[i] and cell[i+1]
    for i in range(0,iMax-1):
        # at cell[i]
        R,RInv,Gam,GamAbs = Apm(i) # matrixes
        Ap = np.dot((np.dot(R,Gam+GamAbs)),RInv)
        # at cell[i+1]
        R,RInv,Gam,GamAbs = Apm(i+1)
        Am = np.dot((np.dot(R,Gam-GamAbs)),RInv)

        fluxPlus[i] = 0.5*(np.dot(Ap,Qc[i]) + np.dot(Am,Qc[i+1]))

def Apm(n):
    H = (Qc[n][2] + qf[n][2])/Qc[n][0] # enthalpy
    u = qf[n][0]
    c = np.sqrt((gamma-1)*(H-0.5*u**2))
    b2 = (gamma-1)/c**2
    b1 = 0.5*b2*u**2

    R = np.array([[1.,    1.,       1.],\
                  [u-c,   u,        u+c],\
                  [H-u*c, 0.5*u**2, H+u*c]])
    RInv = np.array([[0.5*(b1+u/c), 0.5*(-b2*u-1/c), 0.5*b2],\
                     [1-b1,         b2*u,            -b2],\
                     [0.5*(b1-u/c), 0.5*(-b2*u+1/c), 0.5*b2]])
    Gam = np.array([[u-c, 0., 0.],\
                    [0.,  u,  0.],\
                    [0.,  0., u+c]])
    GamAbs = np.array([[abs(u-c), 0.,     0.],\
                       [0.,       abs(u), 0.],\
                       [0.,       0.,     abs(u+c)]])
    return R,RInv,Gam,GamAbs

def qf2Qc(qf): # qf -> Qc
    lo_Qc = np.zeros((iMax,3))
    for i in range(iMax):
        lo_Qc[i][0] = qf[i][1]
        lo_Qc[i][1] = qf[i][1]*qf[i][0]
        lo_Qc[i][2] = (qf[i][2]/(gamma-1) + 0.5*qf[i][1]*(qf[i][0]**2))
    return lo_Qc

def Qc2qf(Qc):
    lo_qf = np.zeros((iMax,3))
    for i in range(iMax):
        lo_qf[i][0] = Qc[i][1]/Qc[i][0]
        lo_qf[i][1] = Qc[i][0]
        lo_qf[i][2] = (gamma-1)*(Qc[i][2]-0.5*Qc[i][0]*((Qc[i][1]/Qc[i][0])**2))
    return lo_qf

def output(dir,file,q,x):
    outData = ["x    u     rho     p"]
    for i in range(len(q)):
        outData.append(str(x[i])+" "+str(q[i][0])+" "+str(q[i][1])+" "+str(q[i][2]))
    outData = "\n".join(outData)

    with open(dir+"/"+file,"wt") as f:
        f.write(outData)

def mkdir():
    try:
        os.mkdir("output")
    except:
        pass

mkdir()
initialCondition()
for i in range(nStep):
    print(i)
    calQ()
    qf = Qc2qf(Qc)

    output("output","time"+"{:0=4}".format(int(i))+"d-3",qf,x)
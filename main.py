
#Gauss seidal solution process for simple power network

#Written By: Amila Nuwan Pathirana, Karasin Pthirannnahalage

# Date: 2014-10-17 modified on 2019-10-20

# import numpy library for complex calculations
import numpy as np

#Function to convert phasor to rectangular form
#Limit the answer to 4 decimal places

decimalplaces=5

def P2R(A, phi):
    R=A*np.cos(phi)
    R = np.around(R, decimals=decimalplaces)
    I= A*np.sin(phi)*1j
    I = np.around(I, decimals=decimalplaces)
    Rect=R+I
    return Rect

#Function to convert rectangular to phasor  form
#Limit the answer to 4 decimal places
def R2P(R):
    r = R.real
    i = R.imag
    Vm = np.sqrt(r ** 2 + i ** 2)
    Vm=np.around(Vm,decimals=decimalplaces)
    Va = np.arctan(i / r)
    Va = np.around(Va, decimals=decimalplaces)
    return Vm,Va

# Y Matrix
Y=[[-7.5j,2.5j,5j],
   [2.5j,-6.5j,4j],
   [5j,4j,-9j]]

#Bus Index
    #0- Slack
    #1-PQ
    #2-PV

Bus=np.array([0,1,2])  # First bus is a slack, 2nd Bus is a PQ bus and 3rd Bus is a PV bus

V = np.ones(len(Bus),dtype=complex)
S = np.zeros(len(Bus),dtype=complex)

# Voltages at each bus, I f dont have any voltage magnitudes or angles we assume a reasonable
# value for those variables
# Unknown voltage magnitudes assume 1pu
# Unknown voltage angles assume 0

V[0]=P2R(1.04,0) # here we know both values
V[1]=P2R(1.0,0) # here we dont know both magnitude and angle. So we assume angle to be zero and
# magnitude to be 1
V[2]=P2R(1.005,0) # here we know the magnitude but don't know the angle. we assume the angle to
# be zero

S[1]= -1-.8j  # P Q injection "INTO" the bus
S[2]=1+0j  #We know the power injection because it is an generator bus. Q we dont know so we put
# Zero for now

# Round the S,V values to 4 decimal places
S=np.around(S,decimals=decimalplaces)
V=np.around(V,decimals=decimalplaces)

max_iterations=10

print("\nGauss seidal solution process for simple power network. Written By: Amila Pathirana")
print("\n////////////////////////////////////////////////////////////")
print("Our final objective is to find all the voltage magnitudes and angles at each and every\n"
      "bus. For slack bus we know the voltage magnitude and angle. for PQ bus we directly \n "
      "calculate V[i] ,for PV buses, first we calculate Q[i] first and then calculate V[i] using\n "
      "the calculated Q[i].")
print("////////////////////////////////////////////////////////////\n")
for x in range(max_iterations):

    for i in range(len(Bus)):
        if Bus[i]==0:
            V[i]=P2R(1.04,0)

        if Bus[i]==1:
            Sum = np.zeros(1,dtype=complex)
            for j in range(len(Bus)):
                Sum+=Y[i][j]*V[j]

            Sum=Sum-Y[i][i]*V[i]
            V[i]=(1/Y[i][i])*((np.conjugate(S[i])/np.conjugate(V[i]))-Sum)

        if Bus[i]==2:
            Vt=V[i]
            Sum2 = np.zeros(1,dtype=complex)
            for j in range(len(Bus)):
                Sum2+=Y[i][j]*V[j]

            Sum4=Sum2*np.conjugate(-1*V[i])

            P=S[i].real
            Q=Sum4.imag
            S[i]=P+Q*(1j)  # keep the P as it is, Change only the Q . Remember this is a PV bus

            Sum3 = np.zeros(1, dtype=complex)
            for j in range(len(Bus)):
                Sum3 += Y[i][j] * V[j]

            Sum3=Sum3-Y[i][i]*V[i]
            # Since we know Q now, we can calculate V[i] as we did for PQ bus
            V[i]=(1/Y[i][i])*((np.conjugate(S[i])/np.conjugate(V[i]))-Sum3)

            Vr=Vt.real
            Vi=Vt.imag
            Vm=np.sqrt(Vr**2+Vi**2)
            Vr=V[i].real
            Vi=V[i].imag
            Va=np.arctan(Vi/Vr)
            V[i] = P2R(Vm, Va) # we keep the voltage magnitude same as the previous step. THIS IS
            # A PV bus.we only change the angle of the voltage

        V[i] = np.around(V[i],decimals=decimalplaces) # Limit decimal places of calculated V

    Vol, Ang = R2P(V)
    Vol = np.around(Vol, decimals=decimalplaces)
    Ang = np.around(Ang * 180 / np.pi, decimals=decimalplaces)  # angle is converted to degrees

    print("Voltage Magnitudes After Iteration ", x + 1, ":", Vol)
    print("Voltage Angles (Deg) After Iteration ", x + 1, ":", Ang)
    print("////////////////////////////////////////////////////////////\n")
    print(" :) ")







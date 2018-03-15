# if you dont have the XRAYLIB package 
# download it here (https://github.com/tschoonj/xraylib/wiki/Installation-instructions)
import xraylib 
import numpy as np


##################################################################
########## DEFINE FUNCTIONS TO DO ANNOYING STUFF##################
##################################################################
def Density(layerdict):# send a string of the compound of interest
    Material = layerdict['Name']
    if Material == 'ZnO':
            return 5.6 #g/cm3
    elif Material == 'CIGS':
            return 5.75 #g/cm3
    elif Material == 'ITO':
            return 7.14 #g/cm3
    elif Material == 'CdS':
            return 4.826 #g/cm3
    elif Material == 'Kapton': # http://physics.nist.gov/cgi-bin/Star/compos.pl?matno=179
            return 1.42 #g/cm3 
    elif Material == 'SiN':
            return 3.44 #g/cm3
    if Material == 'Mo':
            return 10.2 #g/cm3

def GetMaterialMu(E, data): # send in  the photon energy and the dictionary holding the layer information
    Ele = data['Element']
    Mol = data['MolFrac']
    t = 0
    for i in range(len(Ele)):
            t += xraylib.AtomicWeight(xraylib.SymbolToAtomicNumber(Ele[i]))*Mol[i]
    mu=0
    for i in range(len(Ele)):
            mu+= (xraylib.CS_Total(xraylib.SymbolToAtomicNumber(Ele[i]),E) * 
                        xraylib.AtomicWeight(xraylib.SymbolToAtomicNumber(Ele[i]))*Mol[i]/t)
    return mu # total attenuataion w/ coherent scattering in cm2/g

def GetFluorescenceEnergy(Element,Beam): # send in the element and the beam energy to get the Excited Fluorescence Energy 
 #this will return the highest energy fluorescence photon capable of being excited by the beam
    Z = xraylib.SymbolToAtomicNumber(Element)
    F = xraylib.LineEnergy(Z,xraylib.KA1_LINE)
    if xraylib.EdgeEnergy(Z,xraylib.K_SHELL) > Beam:
            F = xraylib.LineEnergy(Z,xraylib.LA1_LINE)
            if xraylib.EdgeEnergy(Z,xraylib.L1_SHELL) > Beam:
                    F = xraylib.LineEnergy(Z,xraylib.LB1_LINE)
                    if xraylib.EdgeEnergy(Z,xraylib.L2_SHELL) > Beam:
                            F = xraylib.LineEnergy(Z,xraylib.LB1_LINE)
                            if xraylib.EdgeEnergy(Z,xraylib.L3_SHELL) > Beam:
                                    F = xraylib.LineEnergy(Z,xraylib.LG1_LINE)
                                    if xraylib.EdgeEnergy(Z,xraylib.M1_SHELL) > Beam:
                                            F = xraylib.LineEnergy(Z,xraylib.MA1_LINE)
    return F

##############################################################
########### USER INPUT VALUES AND I/IO CALCULATION ###########
##############################################################


#define each of the layers that are in the stack, include a: 
# list of elements, 
# the stoicheometry of each element, 
# the thickness of the layer, (microns)
# the layer name (used to find density)
ZNO = {'Element':['Zn','O'],'MolFrac':[1,1],'Thick':0.2,'Name':'ZnO'}
CDS = {'Element':['Cd','S'],'MolFrac':[1,1],'Thick':0.05,'Name':'CdS'}
CIGS = {'Element':['Cu','In','Ga','Se'],'MolFrac':[0.8,0.7,0.3,2],'Thick':1.8,'Name':'CIGS'}
MO = {'Element':['Mo'],'MolFrac':[1],'Thick':0.7,'Name':'Mo'}

#COMBINE THE LAYERS FROM ABOVE INTO A LIST (TOP LAYER FIRST, BOTTOM LAYER LAST)
STACK = [ZNO,CDS,CIGS,MO]

#here you define the measurement setup
Beam_Theta = 75 #degrees
Detector_Theta = 15 #degrees
Beam_Energy = 10.5 #keV

# here you define the element you want the attenuation correction for
EOI = 'Cu'

iio = 1
um_to_cm = 10**-4


for layer1 in STACK:
    if EOI in layer1['Element']:
        t = layer1['Thick']*um_to_cm
        dt = 0.01*um_to_cm # 10 nm stepsizes
        steps = int(t/dt);
        first_integral = [None]*steps
        path_in = np.zeros((steps,1))
        path_out = np.zeros((steps,1))
        for N in range(steps):
            for idx in range(N+1):
                path_in[idx] = -Density(layer1)*GetMaterialMu(Beam_Energy,layer1)*dt/np.sin(Beam_Theta*np.pi/180)
                path_out[idx] = -Density(layer1)*GetMaterialMu(GetFluorescenceEnergy(EOI,Beam_Energy),layer1)*dt/np.sin(Detector_Theta*np.pi/180)
            first_integral[N] = np.exp(np.sum(path_in+path_out))
        iio = iio*np.sum(first_integral)/steps
        break
    else:
        t = layer1['Thick']*um_to_cm
        path_in =  -Density(layer1) * GetMaterialMu(Beam_Energy,layer1) * t/np.sin(Beam_Theta*np.pi/180)
        path_out = -Density(layer1) * GetMaterialMu(GetFluorescenceEnergy(EOI,Beam_Energy),layer1) * t/np.sin(Detector_Theta*np.pi/180)
        iio = iio*np.exp(path_in + path_out) 
            
print(iio)

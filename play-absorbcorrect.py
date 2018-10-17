import xraylib as xl
import numpy as np

def ele_dens(atomicNumber):
    D = xl.ElementDensity(atomicNumber)
    return D

def ele_weight(atomicNumber):
    W = xl.AtomicWeight(atomicNumber)
    return W

def symToNum(sym):
    S = xl.SymbolToAtomicNumber(str(sym))
    return S

def capXsect(element, energy):
    C = xl.CS_Total(symToNum(element), energy)
    return C

#returns the highest energy fluorescence photon of element capable of being excited by beam
    #str() necessary to make list format of elements compatible with byte read of xraylib
def XRF_line(Element,Beam_Energy):
    Z = xl.SymbolToAtomicNumber(str(Element))
    F = xl.LineEnergy(Z,xl.KA1_LINE)
    if xl.EdgeEnergy(Z,xl.K_SHELL) > Beam_Energy:
            F = xl.LineEnergy(Z,xl.LA1_LINE)
            if xl.EdgeEnergy(Z,xl.L1_SHELL) > Beam_Energy:
                    F = xl.LineEnergy(Z,xl.LB1_LINE)
                    if xl.EdgeEnergy(Z,xl.L2_SHELL) > Beam_Energy:
                            F = xl.LineEnergy(Z,xl.LB1_LINE)
                            if xl.EdgeEnergy(Z,xl.L3_SHELL) > Beam_Energy:
                                    F = xl.LineEnergy(Z,xl.LG1_LINE)
                                    if xl.EdgeEnergy(Z,xl.M1_SHELL) > Beam_Energy:
                                            F = xl.LineEnergy(Z,xl.MA1_LINE)
    return F

#calculate total attenuataion w/ coherent scattering (cm2/g) at beam energy for element in layer
def MatMu(energy, layer):
    layer_elements = layer['Element']
    layer_molFrac = layer['MolFrac']
    layer_ele_index = range(len(layer_elements))
    
    layer_element_mol = 0
    for ele in layer_ele_index:
        layer_element_mol += ele_weight(symToNum(layer_elements[ele]))*layer_molFrac[ele]
        #in chosen layer: takes first element atomic weight, multiplies by mol fraction, 
            #takes second element atomic weight, multiplies by mol fraction, 
            #adds to previous element weight, etc. Final number is total g/mol of elements in layer?
    
    layer_ele_mu = 0
    for ele in layer_ele_index:
        layer_ele_mu += capXsect(layer_elements[ele], beam_energy) #* ele_weight(symToNum(layer_elements[ele])) * layer_molFrac[ele] / layer_element_mol
        #print(layer_ele_mu)
    return layer_ele_mu

def getSublayers(layer):
    sublayers= []
    for layer in layers:
        T = layer['Thick'] * (1/10000)          #convert um to cm (1cm/ 10000um)
        sublayers = int(T/dt)  
        key = 'numSublayers'
        layer.setdefault(key, sublayers)      

#iio = (1/T) * sum_T(   exp   (sum_N( (-mu_in(t_i) / cos(theta_in)) * dt + ( -mu_out(t_i) / cos(theta_out)) *dt)))
SNO2 = {'Element':['Sn','O'],'MolFrac':[1,2],'Thick':0.6,'LDensity': 6.85, 'Name': 'SnO2'}
CDS = {'Element':['Cd','S'],'MolFrac':[1,1],'Thick':0.08,'LDensity': 4.82, 'Name': 'CdS'}
CDTE = {'Element':['Cd','Te'],'MolFrac':[1,1],'Thick':5.0,'LDensity': 5.85, 'Name': 'CdTe'}
CU = {'Element':['Cu'],'MolFrac':[1],'Thick':0.01,'LDensity': 8.96, 'Name': 'Cu'}
ZNTE = {'Element':['Zn','Te'],'MolFrac':[1,1],'Thick':0.375,'LDensity': 6.34, 'Name': 'ZnTe'}
MO = {'Element':['Mo'],'MolFrac':[1],'Thick':0.5,'LDensity': 10.2, 'Name': 'Mo'}

#COMBINE THE LAYERS FROM ABOVE INTO A LIST (TOP LAYER FIRST, BOTTOM LAYER LAST)
layers = [MO, ZNTE, CU, CDTE, CDS, SNO2]

beam_energy = 12.7  #keV
beam_theta = 90     #degrees
beam_geometry = np.sin(beam_theta*np.pi/180)
detect_theta = 47                   #degrees
detect_gemoetry = np.sin(detect_theta*np.pi/180)

EOIs = ['Cd', 'Te', 'Cu']
#EOIs = input('Enter elements you wish to correct for in a list and as strings: ')
iio = 1
dt = 10 * (1*10**-3) * (1*10**-4)          #convert 10nm to cm: 10nm * (1um/1000nm) * (1cm/10000um)

for layer in layers:
    print()
    for EOI in EOIs:
        if EOI in layer['Element']:
            getSublayers(layer)
            integral = [None]*layer['numSublayers']
            path_in = np.zeros((layer['numSublayers'],1))
            path_out = np.zeros((layer['numSublayers'],1))
            for sublayer in range(layer['numSublayers']):
                for dx in range(sublayer+1):
                    path_in[dx] = -layer['LDensity'] * MatMu(beam_energy,layer) * dt / beam_geometry
                    path_out[dx] = -layer['LDensity'] * MatMu(XRF_line(EOI,beam_energy), layer) * dt / detect_gemoetry
                integral[sublayer] = np.exp(np.sum(path_in+path_out))
            iio = iio * np.sum(integral)/layer['numSublayers']
            print(EOI, layer['Name'], iio)

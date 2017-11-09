def ABSORB(Beam_Theta,Detector_Theta,Beam_Energy,t):

	import xraylib 
	import numpy as np
	from time import time

	xraylib.XRayInit()
	xraylib.SetErrorMessages(0)

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

	def Density(Material):# send a string of the compound of interest
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
			
	def GetLayerInfo(Layer): #send in a string to get typical layer thickness and dictionary of composition
			um_to_cm = 10**-4
			
			if Layer == 'ZnO':
					mat = {'Element':['Zn','O'],'MolFrac':[1,1]}
					t = 0.2*um_to_cm
					return mat,t
			elif Layer == 'CdS':
					mat = {'Element':['Cd','S'],'MolFrac':[1,1]}
					t = 0.05*um_to_cm
					return mat,t
			elif Layer == 'Kapton':
					mat = {'Element':['H','C','N','O'],'MolFrac':[0.026362,0.691133,0.073270,0.209235]} # http://physics.nist.gov/cgi-bin/Star/compos.pl?matno=179
					t = 26.6*um_to_cm #measured using profilometer
					return mat, t
			elif Layer == 'ITO':
					mat = {'Element':['In','Sn','O'],'MolFrac':[1.8,0.1,2.9]} #90% In2O3 #10% SnO2
					t = 0.15*um_to_cm
					return mat,t
			elif Layer == 'Mo':
					mat = {'Element':['Mo'],'MolFrac':[1]}
					t = 0.7*um_to_cm
					return mat,t

			
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

	def GetIIO(Layer,Energy):
			ROI,t = GetLayerInfo(Layer)
			return np.exp(-Density(Layer)*GetMaterialMu(Energy,ROI)*t)
	
	#conversion factor
	um_to_cm = 10**-4	
	t = t*um_to_cm

	##Set incident Beam Energy and Detector Geometry
	# Beam_Theta = 90 #degrees
	# Detector_Theta = 47 #degrees
	# Beam_Energy = 10.5 #keV
	# Get Layor of interest information
	L = 'CIGS'
	ROI = {'Element':['Cu','In','Ga','Se'],'MolFrac':[0.8,0.8,0.2,2]}
	Elem = ROI['Element']

	# define sublayers thickness and adjust based on measurement geometry
	dt = 0.01*um_to_cm # 10 nm stepsizes
	steps = int(t/dt);
	T = np.ones((steps,1))*dt 
	beam_path = T/np.sin(Beam_Theta*np.pi/180)
	fluor_path = T/np.sin(Detector_Theta*np.pi/180)

	# initialize variables to hold correction factors
	iio = [None]*steps
	factors = [None]*len(Elem)
	print 'For a film thickness of ', t/um_to_cm, 'microns:'
	#loop over sublayers for self attenuation and top layer attenuation
	ti = time()
	for ind,Z in enumerate(Elem):
			for N in range(steps):
					beam_in = -Density(L)*GetMaterialMu(Beam_Energy,ROI)*beam_path[0:N]
					beam_out = -Density(L)*GetMaterialMu(GetFluorescenceEnergy(Z,Beam_Energy),ROI)*fluor_path[0:N]
					iio[N] = np.exp(np.sum(beam_in+beam_out))
			factors[ind] = np.sum(iio)/N * GetIIO('Kapton',Beam_Energy) * GetIIO('Kapton',GetFluorescenceEnergy(Z,Beam_Energy))
			print 'The absorption of ', Z, 'in', L,'at beam energy', Beam_Energy,'is', round(factors[ind]*100,2)
	print 'Calculation Time = ', round(time()-ti,2),'s for ', steps, 'iterations on ', ind+1,' elements'
	return factors
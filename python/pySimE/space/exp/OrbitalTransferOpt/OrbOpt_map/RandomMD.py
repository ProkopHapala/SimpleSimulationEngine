

from pylab import *

def MCMD_Run( fitnesFunc, XBest, stepSize, accuracy, stopAfter, missRatio, maxIter, GenHistory,  wStep=0.7, fBias = 3.0,  kBias0=0.5, kBiasf = 0.95 ):
	fitness = 100000000
	fitness = fitnesFunc(Gen)
	print 
	print " ========= MC-(quasi)MD Optimization ================= "
	print "stopAfter, missRatio:  ", stopAfter, missRatio
	print " Initial === ", Gen, fitness
	fromMajor = 0
	fitnessMajor = fitness
	kBias = 0.0
	V = zeros(len(Gen))
	dx = zeros(len(Gen))
	fitness = 
	for i in range(maxIter):
		n = len(Gen)
		dX    =  V +  (rand(n)-0.5)*stepSize
		dXdX  = dot(dXdX)
		X     =  XBest + (rand(n)-0.5)*stepSize 
		E = fitnesFunc(XNew)
		if(E<EBest):
			print " fitness: ",fitness," step ",stepSize," kBias ", kBias, " i: ", i," fromMajor",fromMajor
			V     += dX
			XBest  = X
			EBest  = E
			if abs(EMajor-EBest)>accuracy:
				fitnessMajor = fitness
				fromMajor    = 0
		else:
			F = (fittness-fitnessBest) * dX / dXdX       #  dE * dX / |dX|^2
			V *= kFrict    # Friction
			V -= F/mass    # Slope force
			V -= fireCoef * ( V - F*sqrt(dot(v,v)/dot(F,F)) )  # Fire MD algorithm
	if fromMajor>stopAfter:
		print " Not able to improve => Exiting .... " 
		break
	fromMajor += 1
	return XBest

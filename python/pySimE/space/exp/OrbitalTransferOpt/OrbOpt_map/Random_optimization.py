



from pylab import *

def TryNew( fitnesFunc, Gen, fitnessBest, stepSize ):
    hit = False
    n = len(Gen)
    GenNew     =  Gen + (rand(n)-0.5)*stepSize 
    fitnessNew = fitnesFunc(GenNew)
    #print GenNew, fitnessNew
    if(fitnessNew < fitnessBest ):
        hit = True
        Gen    = GenNew
        fitnessBest = fitnessNew
    return Gen, fitnessBest,hit

def TryNew2( fitnesFunc, Gen, fitnessBest, stepSize ):
    hit = False
    n = len(Gen)
    GenNew     =  Gen + (rand(n)-0.5)*stepSize 
    fitnessNew = fitnesFunc(GenNew)
    #print GenNew, fitnessNew
    if(fitnessNew < fitnessBest ):
        hit = True
        Gen    = GenNew
        fitnessBest = fitnessNew
    return Gen, fitnessBest,hit

def MC_Run( fitnesFunc, Gen, stepSize, accuracy, stopAfter, missRatio, maxIter, GenHistory ):
    fitness = 100000000
    fitness = fitnesFunc(Gen)
    print 
    print " ========= MC Optimization ================= "
    print "stopAfter, missRatio:  ", stopAfter, missRatio
    print " Initial === ", Gen, fitness
    badInRow  = 0
    fromMajor = 0
    fitnessMajor = fitness
    for i in range(maxIter):
        Gen,fitness,hit = TryNew( fitnesFunc, Gen, fitness, stepSize )
        if(hit):
            #print " fitness: ",fitnessBest, " stepSize: ", stepSize 
            badInRow = 0
            stepSize *= 1.1
            print " fitness: ",fitness," step ",stepSize, " i: ", i, " bad ",badInRow," fromMajor",fromMajor
            GenHistory.append(Gen)
            if abs(fitnessMajor-fitness)>accuracy:
                fitnessMajor = fitness
                fromMajor    = 0
                #print " === restart fromMajor "
        if badInRow>missRatio:
            stepSize *= 0.75
            badInRow = 0
            print " stepSize down to ", stepSize
        if fromMajor>stopAfter:
            print " Not able to improve => Exiting .... " 
            break
        badInRow  += 1
        fromMajor += 1
    return Gen
    
def MCBias_Run( fitnesFunc, Gen, stepSize, accuracy, stopAfter, missRatio, maxIter, GenHistory ):
    fitness = 100000000
    fitness = fitnesFunc(Gen)
    print 
    print " ========= MC Bias Optimization ================= "
    print "stopAfter, missRatio:  ", stopAfter, missRatio
    print " Initial === ", Gen, fitness
    badInRow  = 0
    fromMajor = 0
    fitnessMajor = fitness
    kBias = 0.0
    dGenNorm = 0.0
    dGen = zeros(len(Gen))
    for i in range(maxIter):
        GenNew,fitness,hit = TryNew( fitnesFunc, Gen + kBias*dGen, fitness, (stepSize+kBias*dGenNorm)*(1.0-kBias) )
        if(hit):
            dGen = (GenNew - Gen) * 1.2
            dGenNorm = dot(dGen,dGen)
            #dGen /= dot(dGen,dGen)
            Gen = GenNew
            #print " fitness: ",fitnessBest, " stepSize: ", stepSize 
            badInRow = 0
            stepSize *= 1.1
            kBias = 0.9
            print " fitness: ",fitness," step ",stepSize, " i: ", i, " bad ",badInRow," fromMajor",fromMajor
            GenHistory.append(Gen)
            if abs(fitnessMajor-fitness)>accuracy:
                fitnessMajor = fitness
                fromMajor    = 0
                #print " === restart fromMajor "
        else:
            kBias *= 0.75
        if badInRow>missRatio:
            stepSize *= 0.75
            badInRow = 0
            print " stepSize down to ", stepSize
        if fromMajor>stopAfter:
            print " Not able to improve => Exiting .... " 
            break
        badInRow  += 1
        fromMajor += 1
    return Gen

def MCBias2_Run( fitnesFunc, Gen, stepSize, accuracy, stopAfter, missRatio, maxIter, GenHistory,  wStep=0.7, fBias = 3.0,  kBias0=0.5, kBiasf = 0.95 ):
    fitness = 100000000
    fitness = fitnesFunc(Gen)
    print 
    print " ========= MC Bias Optimization ================= "
    print "stopAfter, missRatio:  ", stopAfter, missRatio
    print " Initial === ", Gen, fitness
    fromMajor = 0
    fitnessMajor = fitness
    kBias = 0.0
    dGen = zeros(len(Gen))
    for i in range(maxIter):
        GenNew,fitness,hit = TryNew( fitnesFunc, Gen + fBias*kBias*dGen, fitness, stepSize*(1.0-kBias) )
        if(hit):
            print " fitness: ",fitness," step ",stepSize," kBias ", kBias, " i: ", i," fromMajor",fromMajor
            dGen     = GenNew - Gen
            stepSize = wStep*stepSize + (1.0-wStep)*sqrt(dot(dGen,dGen))
            kBias    = kBias0
            Gen      = GenNew
            GenHistory.append(Gen)
            if abs(fitnessMajor-fitness)>accuracy:
                fitnessMajor = fitness
                fromMajor    = 0
                #print " === restart fromMajor "
        else:
            stepSize *= kBiasf
            kBias    *= kBiasf**4
            #print " ------- stepSize ",stepSize," kBias  ", kBias
        if fromMajor>stopAfter:
            print " Not able to improve => Exiting .... " 
            break
        fromMajor += 1
    return Gen
				
def MCBias2b_Run( fitnesFunc, Gen, stepSize, accuracy, stopAfter, missRatio, maxIter, GenHistory,  kElong=1.3, kContr=0.95, fBias = 3.0,  kBias0=0.5, kBiasf = 0.95 ):
    fitness = 100000000
    fitness = fitnesFunc(Gen)
    print 
    print " ========= MC Bias Optimization ================= "
    print "stopAfter, missRatio:  ", stopAfter, missRatio
    print " Initial === ", Gen, fitness
    fromMajor = 0
    fitnessMajor = fitness
    kBias = 0.0
    dGen = zeros(len(Gen))
    for i in range(maxIter):
        GenNew,fitness,hit = TryNew( fitnesFunc, Gen + dGen, fitness, stepSize*(1.0-kBias) )
        if(hit):
            print " fitness: ",fitness," step ",stepSize," kBias ", kBias, " i: ", i," fromMajor",fromMajor
            dGen      = GenNew - Gen
            stepSize *= kElong
            kBias    *= kBias0
            Gen       = GenNew
            GenHistory.append(Gen)
            if abs(fitnessMajor-fitness)>accuracy:
                fitnessMajor = fitness
                fromMajor    = 0
                #print " === restart fromMajor "
        else:
            stepSize *= kBiasf
            kBias    *= kBiasf**4
            #print " ------- stepSize ",stepSize," kBias  ", kBias
        if fromMajor>stopAfter:
            print " Not able to improve => Exiting .... " 
            break
        fromMajor += 1
    return Gen

def MCBias3_Run( fitnesFunc, Gen, stepSize, accuracy, stopAfter, missRatio, maxIter, GenHistory ):
    fitness = 100000000
    fitness = fitnesFunc(Gen)
    print 
    print " ========= MC Bias Optimization ================= "
    print "stopAfter, missRatio:  ", stopAfter, missRatio
    print " Initial === ", Gen, fitness
    fromMajor = 0
    fitnessMajor = fitness
    kBias = 0.0
    dGen = zeros(len(Gen))
    for i in range(maxIter):
        GenRnd     =  rand(n)-0.5
        GenRnd     =  GenRnd / sqrt( dot(GenRnd,GenRnd) )
        GenNew     =  Gen + GenRnd*stepSize + dGen
    	fitnessNew = fitnesFunc(GenNew)
        GenNew,fitness,hit = TryNew( fitnesFunc, Gen + 1.6*kBias*dGen, fitness, stepSize*(1.0-kBias) )
        if(hit):
            dGen     = GenNew - Gen
            stepSize = sqrt(dot(dGen,dGen))
            kBias = 0.9
            Gen      = GenNew
            #print " fitness: ",fitnessBest, " stepSize: ", stepSize 
            print " fitness: ",fitness," step ",stepSize, " i: ", i," fromMajor",fromMajor
            GenHistory.append(Gen)
            if abs(fitnessMajor-fitness)>accuracy:
                fitnessMajor = fitness
                fromMajor    = 0
                #print " === restart fromMajor "
        else:
            kBias *= 0.9
        if fromMajor>stopAfter:
            print " Not able to improve => Exiting .... " 
            break
        fromMajor += 1
    return Gen

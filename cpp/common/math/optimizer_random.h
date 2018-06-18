
typedef double (*FitnessFunction)( double * xs );

class OptimizerRandom{
	public:
		int n     = 0;
		int nfail = 0;
		double * xBest = NULL;
		double * xNew  = NULL;
		double * dxMax = NULL;
		double bestFitness = 0.0f;
		FitnessFunction fitnessFunction;

	virtual void restart(){
        bestFitness = fitnessFunction( xBest );
        nfail = 0;
	}

	OptimizerRandom(int n_, double* xBest_, double* dxMax_, FitnessFunction fitnessFunction_  ){
		n = n_;
		dxMax = dxMax_;
		xBest = new double [n];	for(int i=0; i<n; i++){  xBest[i] = xBest_[i]; }
		xNew  = new double [n];
		fitnessFunction = fitnessFunction_;
		restart();
	};

	virtual double step(){
		for(int i=0; i<n; i++){    xNew[i] = xBest[i] + dxMax[i] * ( 2.0*randf() - 1.0 );  }
		double dfitness  = fitnessFunction( xNew ) - bestFitness;
		if( dfitness > 0 ){
			for(int i=0; i<n; i++){ xBest[i] = xNew[i]; }
			bestFitness += dfitness;
			nfail=0;
		}else{
            nfail++;
		}
		return dfitness;
	};

};


class OptimizerRandom_2 : public OptimizerRandom {
	public:
		double decay      = 0.0d;
		double c_moment   = 0.0d;
		double * momentum = NULL;

	OptimizerRandom_2  (int n_, double* xBest_, double* dxMax_, double decay_, double c_moment_, FitnessFunction fitnessFunction_  )
	: OptimizerRandom  ( n_, xBest_, dxMax_, fitnessFunction_  ) {
		decay = decay_;
		c_moment = c_moment_;
		momentum  = new double [n]; for(int i=0; i<n; i++){  momentum[i] = 0; }
	};

	void cleanMomentum(){ for(int i=0; i<n; i++){ momentum[i]=0; } }

    virtual void restart(){
        OptimizerRandom::restart();
        cleanMomentum();
	}

	virtual double step(){
		//for(int i=0; i<n; i++){   xNew[i] = xBest[i] + dxMax[i] * ( 2.0*randf() - 1.0 ) + momentum[i];  }
		for(int i=0; i<n; i++){   xNew[i] = xBest[i] + dxMax[i] * ( 2.0*randf() - 1.0 ) + momentum[i]*randf();  }
		double dfitness  = fitnessFunction( xNew ) - bestFitness;
		if( dfitness > 0 ){
			for(int i=0; i<n; i++){
				momentum[i] += c_moment * ( xNew[i] - xBest[i] );
				xBest[i]     = xNew[i];
			}
			bestFitness += dfitness;
		} else {
            nfail++;
			for(int i=0; i<n; i++){ momentum[i] *= decay; }
		}
		return dfitness;
	};

};




class OptimizerRandom_3 : public OptimizerRandom {
	public:
		double decay      = 0.0d;
		double c_moment   = 0.0d;
		double * xBestNew = NULL;
		double bestNewFitness;

	OptimizerRandom_3  (int n_, double* xBest_, double* dxMax_, double decay_, double c_moment_, FitnessFunction fitnessFunction_  )
	: OptimizerRandom  ( n_, xBest_, dxMax_, fitnessFunction_  ) {
		decay = decay_;
		c_moment = c_moment_;
		xBestNew  = new double [n]; for(int i=0; i<n; i++){  xBestNew[i] = xBest[i]; }
	};

	void success_step(){
		for(int i=0; i<n; i++){
			xBestNew[i]  = xNew[i] + c_moment * ( xNew[i] - xBest[i] );
			xBest[i]     = xNew[i];
		}
		bestNewFitness = fitnessFunction( xNew );
		if( bestNewFitness > bestFitness ){
			bestFitness  = bestNewFitness;
			success_step();
		}
	}

	virtual double step(){
		for(int i=0; i<n; i++){   xNew[i] = xBestNew[i] + dxMax[i] * ( 2.0*randf() - 1.0 );  }
		double  fitness = fitnessFunction( xNew );
 		double dfitness = fitness - bestFitness;
		if( dfitness > 0 ){
			bestFitness    = fitness;
			success_step();
		} else {
            nfail++;
			if( fitness > bestNewFitness ){
				bestNewFitness = fitness;
				for(int i=0; i<n; i++){ xBestNew[i]    = xNew[i]; }
			}else{
				for(int i=0; i<n; i++){ xBestNew[i] = decay*xBestNew[i] + (1.0d-decay)*xBest[i]; }
			}
		}
		return dfitness;
	};

};

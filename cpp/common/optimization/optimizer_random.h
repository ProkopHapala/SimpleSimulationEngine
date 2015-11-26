
typedef double (*FitnessFunction)( double * xs );


int nEvals = 0;

double fittnessFunc_rosenbrok( double * xs ){
	nEvals++;
	//return -rosenbrok( xs[0], xs[1] );
	//return -sinValey( xs[0], xs[1] );
	//return -cosValey( xs[0], xs[1] );
	return -spiral( xs[0], xs[1] );
}

class OptimizerRandom{
	public:
		int n;
		double * xBest ;
		double * xNew  ;
		double * dxMax ;
		double bestFitness;
		FitnessFunction fitnessFunction;

	OptimizerRandom(int n_, double* xBest_, double* dxMax_, FitnessFunction fitnessFunction_  ){
		n = n_;
		dxMax   = dxMax_;
		xBest = new double [n];	for(int i=0; i<n; i++){  xBest[i] = xBest_[i]; }
		xNew    = new double [n];		
		fitnessFunction = fitnessFunction_;
		bestFitness = fitnessFunction( xBest );
	};

	virtual double step(){
		for(int i=0; i<n; i++){    xNew[i] = xBest[i] + dxMax[i] * ( 2.0*randf() - 1.0 );  }
		double dfitness  = fitnessFunction( xNew ) - bestFitness;
		if( dfitness > 0 ){ 
			for(int i=0; i<n; i++){ xBest[i] = xNew[i]; } 
			bestFitness += dfitness;
		}
		return dfitness;
	};

};


class OptimizerRandom_2 : public OptimizerRandom {
	public:
		double decay;
		double c_moment;
		double * momentum;
		

	OptimizerRandom_2  (int n_, double* xBest_, double* dxMax_, double decay_, double c_moment_, FitnessFunction fitnessFunction_  )
	: OptimizerRandom  ( n_, xBest_, dxMax_, fitnessFunction_  ) {
		decay = decay_;
		c_moment = c_moment_;
		momentum  = new double [n]; for(int i=0; i<n; i++){  momentum[i] = 0; }
	};

	virtual double step(){
		for(int i=0; i<n; i++){   xNew[i] = xBest[i] + dxMax[i] * ( 2.0*randf() - 1.0 ) + momentum[i];  }
		double dfitness  = fitnessFunction( xNew ) - bestFitness;
		if( dfitness > 0 ){ 			
			for(int i=0; i<n; i++){ 
				momentum[i] += c_moment * ( xNew[i] - xBest[i] );
				xBest[i]     = xNew[i]; 
			} 
			bestFitness += dfitness;
		} else {
			for(int i=0; i<n; i++){ momentum[i] *= decay; } 
		}
		return dfitness;
	};

};




class OptimizerRandom_3 : public OptimizerRandom {
	public:
		double decay;
		double c_moment;
		double * xBestNew;
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

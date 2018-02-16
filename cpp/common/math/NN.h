
#ifndef  NN_h
#define  NN_h

//inline dot( int n, double* a, double* b ){ double sum=0; for(int i=0; i<n; i++){ sum+=a[i]*b[i]; }; }

class NN_layer {  public:

    public:
    int nIn,nOut;
    int af=1;                                               // activation function number
    double *xs;                                             // node values of previous layer = input nodes
    double *ys;                                             // node values of this layer
    double *dfs;                                            // derivatives of nonlinar function at Wx
    //double *dys;                                            // full derivative of this layer with respect to layer 0
    //double *dxs;                                            // full derivative of previous layer with respect to layer 0
    double *bias;                                           // bias values
    double *weight;                                         // weights [nIn,nOut]

    //inline double getLinearNode(double* wi){

    inline double doti(int i, double * x ){
        double  Wx = 0.0;
        double* wi = weight + i*nIn;
        for(int j=0; j<nIn; j++){ Wx+=wi[j]*x[j]; };
        return Wx;
    }

    /*
    inline double doti(int i,){
        double Wx=0.0;
        double wi = weight + i*nIn;
        for(int j=0; j<nIn; j++){ Wx+=wi[j]*x[j]; };
        return Wx + bias[i];
    }
    */
    /*
    inline double getWxi(int i){
        double Wx=0.0;
        double wi = weight + i*nIn;
        for(int j=0; j<nIn; j++){ Wx+=wi[j]*x[j]; };
        return Wx + bias[i];
    }
    */

    /*
    // cannot do this way - must recalculate for each input weight
    inline void  getDWxi(int i, double& Wx, double& dWx ){
        Wx=0.0; dWx=0.0;
        double wi = weight + i*nIn;
        for(int j=0; j<nIn; j++){
            double wij = wi[j];
            Wx  += wij*x [j];
            dWx += wij*dx[j];
        };
        Wx + bias[i];
    }
    */

    inline double* getWi(int i){ return weight + i*nIn; };

    void allocate( int nIn, int nOut ){
        //dfs    = new double[nOut];
        bias   = new double[nOut];
        weight = new double[nOut*nIn];
    }

    void initRandom( double wmin, double wmax ){
        int ii=0;
        for(int i=0; i<nOut; i++){
            bias[i]   = randf(wmin,wmax);
            for(int j=0; j<nIn; j++){
                weight[ii] = randf(wmin,wmax);
                ii++;
            }
        }
    }

    void eval(){
        double Wx,y;
        switch(af){
            case 0:
                for(int i=0; i<nOut; i++){
                    ///Wx   = getWxi(i);
                    //Wx    = doti(nIn,getWi(i),xs) + bias[i];
                    Wx      = doti(i,xs) + bias[i];
                    ys  [i] = Wx;
                    dfs [i] = 0;
                }
                break;
            case 1:
                for(int i=0; i<nOut; i++){
                    //Wx   = dot(nIn,getWi(i),xs) + bias[i];
                    Wx     = doti(i,xs) + bias[i];
                    y      = tanh(Wx);
                    ys [i] = y;
                    dfs[i] = 1 - y*y;
                }
                break;
        }
    }

    /*
    void eval(){  // with full chain derivatives
        double Wx,dWx,y,df;
        switch(af){
            case 0:
                for(int i=0; i<nOut; i++){
                    getDWxi(i,Wx,dWx);
                    //df    = 1;
                     ys [i] = Wx;
                    //dys [i] = dWx;
                    //dfs[i] = 0;
                }
                break;
            case 1:
                for(int i=0; i<nOut; i++){
                    getDWxi(i,Wx,dWx);
                     y     = tanh(Wx);
                    df     = 1 - y*y;
                     ys[i] = y;
                    //dys[i] = dWx*df;
                    //dfs[i] = 1 - y*y;
                }
                break;
        }
    }
    */

};

class NN{ public:

    int nLayers  =0;
    int nNodeMax =0;
    NN_layer * layers;

    double *tIn,*tOut;  // temporary arrays for calculation of derivatives

    void eval() {
        for(int i=0;i<nLayers;i++){ layers[i].eval(); }
    }

    void eval_dEdG(){
        for(int il=0; il<nLayers; il++) {            // loop over layers
            NN_layer& L = layers[il];
            for(int i=0; i<L.nOut; i++) {   // loop over
                //double dWx = dot( L.nIn, L.getWi(i),tIn );
                double dWx = L.doti(i,tIn);
                tOut[i]    = dWx * L.dfs[i];
            }
            if( il<(nLayers-1) ){ double * tmp=tIn; tIn = tOut; tOut=tmp; } // swap temp arrays
        }
    }

    void init( int n, int* ns, double* xs, double* ys ){
        nLayers = n;
        layers  = new NN_layer[nLayers];
        layers[0].xs =  xs;
        //layers[0].xs = dxs;
        nNodeMax = ns[0];
        for(int i=0;i<nLayers;i++){
            int nIn  = ns[i  ];
            int nOut = ns[i+1];
            if(nOut>nNodeMax) nNodeMax=nOut;
            layers[i].allocate( nIn, nOut );
            layers[i].initRandom( -1.0,1.0 );
            if(i<(nLayers-1)){
                double*  vs = new double[ nOut ];
                //double* dvs = new double[ n ];
                layers[i  ].ys  =  vs;
                layers[i+1].xs  =  vs;
                //layers[i  ].dys = dvs;
                //layers[i+1].dxs = dvs;
            }
        }
        layers[nLayers-1].ys =  ys;
        //layers[nLayers-1] = dys;

        tIn  = new double[nNodeMax];
        tOut = new double[nNodeMax];
    }

    /*

        for(int i=0; i<nnode[0]; i++) {
            tIn[i] = dfdx[0][i] * L[0]->weight[k][i];
        }
        for(int l=0; l<nlf-1; l++) {            // loop over layers
            for(int i=0; i<nnode[l+1]; i++) {   // loop over
                double sum = 0.0;
                for(int j=0; j<nnode[l]; j++) {
                    sum += L[l+1]->weight[j][i] * tIn[j];
                }
                tOut[i] = sum * dfdx[l+1][i];
                double * tmp; tIn = tOut; tOut=tmp;
                //if(l<nlf-2) tmpinner[l+1][i] = tmpouter[l][i];
            }
        }
        dEdG[k] = tmpouter[nlf-2][0];
    */


    /*
    void PairRuNNer::RuNNer_forces::calc_dEdG(RuNNer_layer **L) {

        double **tmpinner, **tmpouter;

        tmpinner = new double*[nlf-1];
        tmpouter = new double*[nlf-1];
        for(int i=0; i<nlf-1; i++) {
            tmpinner[i] = new double[nnode[i]];
            tmpouter[i] = new double[nnode[i+1]];
        }

        for(int k=0; k<nsym; k++) {
            for(int i=0; i<nnode[0]; i++) {
                tmpinner[0][i] = dfdx[0][i] * L[0]->weight[k][i];
                if(par.normnodes) tmpinner[0][i] /= L[0]->nprevnode;
            }
            for(int l=0; l<nlf-1; l++) {            // loop over layers
                for(int i=0; i<nnode[l+1]; i++) {   // loop over
                    //tmpouter[l][i] = 0.0;
                    double sum = 0.0;
                    for(int j=0; j<nnode[l]; j++) {
                        sum += L[l+1]->weight[j][i] * tmpinner[l][j];
                    }
                    tmpouter[l][i] = sum * dfdx[l+1][i];
                    if(par.normnodes) tmpouter[l][i] /= L[l+1]->nprevnode;
                    if(l<nlf-2) tmpinner[l+1][i]      = tmpouter[l][i];
                }
            }
            dEdG[k] = tmpouter[nlf-2][0];
        }

        for(int i=0; i<nlf-1; i++) {
            delete[] tmpinner[i];
            delete[] tmpouter[i];
        }
        delete[] tmpinner;
        delete[] tmpouter;

        return;

    }
    */
};

#endif




	class AtomTypes{
		public:
		int ntypes;
		int    * Zs;
		char  ** names; 
		double * vdwRs;
		double * vdwEs;
		Uint32 * colors;

		bool loadFromFile( char const* filename ){
			printf(" loading types from: %s \n", filename );
			FILE * pFile;
			pFile = fopen (filename,"r");
			fscanf (pFile, "%i", &ntypes);
			//printf("ntypes %i \n", ntypes );
			Zs       = new    int[ ntypes ];
			vdwRs    = new double[   ntypes ];
			vdwEs    = new double[   ntypes ];
			names    = new char* [   ntypes ];
			colors   = new Uint32[   ntypes ];

			char hexstring[8];
			for (int i=0; i<ntypes; i++){
				names[i] = new char[6];
				fscanf (pFile, " %lf %lf %i %s %s", &vdwRs[i], &vdwEs[i], &Zs[i], names[i], hexstring );
				colors[i] = (Uint32)strtol(hexstring, NULL, 16);
				//printf( "%i %f %f %i %s  %s %i\n", i, vdwRs[i], vdwEs[i], Zs[i], names[i], hexstring, colors[i] );
			}
			fclose (pFile);
			return 0;
		}

		AtomTypes( char const* filename ){  loadFromFile( filename );  };	
	};




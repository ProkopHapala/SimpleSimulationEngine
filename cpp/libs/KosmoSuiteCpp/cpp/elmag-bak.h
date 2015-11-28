

namespace ElMag{

	void biot_savart_element( const Vec3d& R, const Vec3d& dI, Vec3d& B ){
		Vec3d dB;
		dB.set_cross( dI, r );
		double r2 = R.norm2();
		B.add_mul( dB, 1e-7 / ( r2 * sqrt(r2) ) ); 
	}

	void field_at( const Vec3d& where, Vec3d& B, int n, Vec3d * ps, const Vec3d& dIs ){
		for (int i=0; i<n; i++){ biot_savart_element( where - ps[i], dIs[i], B ); }
	}

	void sample_filed ( int m, Vec3d * wheres, Vec3d * Bs, int n, Vec3d * ps, const Vec3d& dIs ){
		for (int i=0; i<m; i++){ Vec3d B; B.set(0.0); field_at( wheres[i], B, n, ps, dIs );  Bs[i].set(B); }
	}
	
}


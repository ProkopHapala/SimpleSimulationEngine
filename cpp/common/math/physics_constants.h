
#ifndef physics_constants_h
#define physics_constants_h

// TODO : this should go elsewhere (physical constants or something)
const double const_hbar_SI      = 1.054571817e-34;    ///< [J.s]  #6.582119569e-16 # [eV/s]
const double const_Me_SI        = 9.10938356e-31;     ///< [kg]
const double const_e_SI         = 1.602176620898e-19; ///< [Coulomb]
const double const_eps0_SI      = 8.854187812813e-12; ///< [F.m = Coulomb/(Volt*m)]
const double const_eV_SI        = 1.602176620898e-19; ///< [J]
const double const_Angstroem_SI = 1.0e-10;

//const double const_K_SI   =  const_hbar_SI*const_hbar_SI/const_Me_SI; // this is wrong
const double const_K_SI   =  const_hbar_SI*const_hbar_SI/(2*const_Me_SI);   // this is correct see schroedinger equation : https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#Preliminaries
const double const_El_SI  =  const_e_SI*const_e_SI/(4.*M_PI*const_eps0_SI);
const double const_Ry_SI  = 0.5 * const_El_SI*const_El_SI/const_K_SI;

const double const_Ry_eV  = 13.6056925944;
//const double const_El_eVA = const_El_SI/( const_e_SI*const_Angstroem_SI );
// Derivation of const_K_eVA from https://en.wikipedia.org/wiki/Hartree
//      Eh= 2Ry = me * (e^2/(4pi eps0 hbar ))^2 = (me/(hbar^2)) * (e^2/(4pi eps0))^2 = (2*me/(hbar^2))/2 * const_El^2 = (const_El^2) / (2*const_K)
//      2Eh=4Ry =  const_El^2/ const_K
//      const_K = const_El^2 / ( 4*Ry ) 
//const double const_K_eVA  = (const_El_eVA*const_El_eVA)/(4*const_Ry_eV);   // this is correct
//const double const_K_eVA  = (const_El_eVA*const_El_eVA)/(2*const_Ry_eV); // this is wrong

const double const_El_eVA = const_El_SI/( const_eV_SI*const_Angstroem_SI                    );
const double const_K_eVA  = const_K_SI /( const_eV_SI*const_Angstroem_SI*const_Angstroem_SI );
//const double const_Ke_eVA = const_K_eVA*1.5;

#endif




/*
 * =====================================================================================
 *
 *       Filename:  sirah_radsup.h
 *
 *    Description:  Sirahna program
 *
 *        Version:  1.0
 *        Created:  02/23/2018 17:31:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wildan Abdussalam (), edun_as@icloud.com
 *   Organization:  Zazi Insitute
 *
 * =====================================================================================
 */

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "mt19937ar.h"
#include "sys/time.h"
#include "gsl/gsl_odeiv2.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_matrix.h"

//#include "histogram.h"

struct ODE_params {
    
    double Delta0;
};

struct timeval tv;
time_t nowtime;
struct tm * nowtm;
char tmbuf[64], buf[64];

/*=====================Indices=========================*/
unsigned atom_i;                  //atom ith 
unsigned atom_j;                  //atom jth
unsigned sigma_down;             //operator sigma_down
unsigned sigma_up;               //operator sigma_up
unsigned Natoms;                 //the number of the atoms
unsigned state;                  //atom's states
unsigned atoms_level;            //atom's state level
unsigned bra;                    //bra
unsigned ket;                    //ket
unsigned Nphotons;               //Photon number
unsigned photon;
unsigned basis;
long unsigned Nstates;           //the number of dot's states
long unsigned Nbasis;            //basis number
unsigned Nrealisation;           //the number of realisation
double cav_area;                 //cavity's area
double atom_density;             //atom's density
double Lx;                       //mesa's length
double Ly;                       //mesa's width
double rad;                      //inter-atom distancess
unsigned realisation;            //realisation
unsigned iter;
unsigned dim;
double hbar;
double gdeph;


/*================== Independent variables ============*/
double r_min;                    //minimum inter-dot distance
double rt;                       //parameter for short-coupling
double n_refr;                   //Refractive index
double lambda;                   //atom's wavelength
double sigma;                    //mismatch standard deviation
double ssigma;
double dt;                       //time_step
double tmax;                     //maximum time propagation
double t;                        //time
unsigned pulselength;            //pulselength of laser
double V;                        //interaction strength
double V0;                       //short-range interaction strength
double delta;                    //laser detuning
double gcav;                     //cavity strength
double G0;                       //G0
double kappa;
double R;
int notimesteps;            //maxtime step
int change;
int cut_off;
double h;
int Niter;
double t1;

/*===================== vector ========================*/
unsigned *excited_num;  //storing the number Rydberg atoms
double *mismatch;       //mismatch energies
double *atom_x;          //x-absis of atom's coordinates
double *atom_y;          //y-ordinat of atom's coordinates
double *y;              //basis

/*=============== dephasing parameters ================*/
double friction;        //meant effect
double Gamma;           //Damping
double noise;           //global noise
double initnoise;       //initialize noise
double dW;              //white noise processes

/*=============== Observables ================*/
double Nryd;            //<n>
double norm;            //normalisation
double Nryd2;           //<n2>
double numerator;       //numerator for g2
double denom1;          //denominator for g2 atom_i
double denom2;          //denominator for g2 atom_j
double prob;            //normalisation

/*=============== file and function ================*/
FILE *outfile;
FILE *outfile1;
char filename[60];
char filenamestate[60];
void init_params_atoms(void);
void init_params_srclock(void);
void init_ground_atom (void);
void init_inverted_atom (void);
void free_alloc_atoms(void);
void init_configuration_atoms(void);
void counting_excited_state_atoms(void);
double rand_gauss(void);
double V_hop (double rij);
double rad_decay(double rij);
void evolution_gsl_nonadaptive_atom (double *detuning);
void evolution_gsl_adaptive_atom (double *detuning);
void write_data_atom(int run);
void checking_rad(void);
void init_symmetry_atom(void);
void see_coupling (void);
void twodrandposition (void);


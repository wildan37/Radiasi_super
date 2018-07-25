/*
 * =====================================================================================
 *
 *       Filename:  mimitian_radsup.c
 *
 *    Description:  Mimitian program
 *
 *        Version:  1.0
 *        Created:  02/23/2018 17:29:32
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wildan Abdussalam (), edun_as@icloud.com
 *   Organization:  Zazi Insitute
 *
 * =====================================================================================
 */

#include "sirah_radsup.h"

void init_params_atoms(void){
    
    atoms_level  = 2;
    Nphotons = 1;
    Nstates = pow(atoms_level, Natoms + Nphotons);      //the number of states
    Nbasis = Nstates * Nstates * 2;                     //the number of basis for DM
    Nrealisation = 1e0;                                 //realisation number
    dim = 1;                                            //Dimensions
    
/* =================================================== */
/* we calculate single emitter momentum
     k0 = nE / (\hbar c).
     in the calculation, we convert to wavelength formula
     k0 = 2pi n / \lambda
     here n = n_refr and \lambda = lambda
*/
    lambda=479.0;                                       //unit (nm)
    n_refr=2.6;                                         //refractive index of CdSe/ZnSe atoms
    gcav = 1;                                           //atom-cavity coupling
    kappa = 0.1;                                        //cavity loss
    R = 0.1;                                            //pumping rate
/* =================================================== */
    ssigma= sigma*1e3/0.658;                            //standard deviation (ns-1) ~ MHz
    atom_density = 1e-3;                                //atom surface density atoms / nm^2
    Lx = sqrt(Natoms/atom_density);                      //V = m / density
    Ly = Lx;                                            //we consider square lattice
    r_min=10.0;                                         //minimum Qds interdistance (nm)
    V0 = 5.0e3/0.658;                                   //(ns-1) ~ MHz
    rt=15.0;                                            // range r_t (nm)
    G0 = 1.0/0.5;                                       //already (ns-1) ~ fast rad decay
    hbar = 6.58211928e-3;                               //planck constant
    gdeph = 0*ssigma;                                  //dephasing
    delta = 0.0;                                        //detuning
    
    ///Free parameters for normal timescale
    Gamma = 1.0;
    friction = 10;
    
    ///Time propagation
    tmax = 2.0e0;
    dt = 1.0e-6;        //second
    notimesteps = tmax/dt;
    printf("%d\n", notimesteps);
    
    //store vector
    excited_num = malloc (Nstates*sizeof(unsigned));        //excited number
    y = malloc (Nbasis*sizeof(double));                     //basis
    mismatch = malloc (Natoms * sizeof(double));             //mismatch energies
    atom_x = (double*)malloc(Natoms *sizeof(double*));        //coordinate_x
    atom_y = (double*)malloc(Natoms *sizeof(double*));        //coordintae_y
        
}

void init_params_srclock (void) {

    atoms_level  = 2;
    Nphotons = 1;
    Nstates = pow(atoms_level, Natoms + Nphotons);      //the number of states
    Nbasis = Nstates * Nstates * 2;                     //the number of basis for DM
    Nrealisation = 1e0;                                 //realisation number
    dim = 1;                                            //Dimensions
    
/* =================================================== */
/* we calculate single emitter momentum
     k0 = nE / (\hbar c).
     in the calculation, we convert to wavelength formula
     k0 = 2pi n / \lambda
     here n = n_refr and \lambda = lambda
*/
///All parameters are scaled by kappa

    lambda=479.0;                                       //unit (nm)
    n_refr=2.6;                                         //refractive index of CdSe/ZnSe atoms
    gcav = 0;                                           //atom-cavity coupling
    kappa = 0.1;                                        //cavity loss
    R = 0.1;                                            //pumping rate
/* =================================================== */
    ssigma= sigma*1e3/0.658;                            //standard deviation (ns-1) ~ MHz
    atom_density = 1e-3;                                //atom surface density atoms / nm^2
    Lx = sqrt(Natoms/atom_density);                      //V = m / density
    Ly = Lx;                                            //we consider square lattice
    r_min=10.0;                                         //r_min - ignore it (nm)
    V0 = 5.0e2/0.658;                                   //(ns-1) ~ MHz
    rt = 15.0;                                          // range r_t (nm)
    G0 = 2.0;                                           //decay constant
    hbar = 6.58211928e-3;                                         //planck constant
    gdeph = 0;                                         //dephasing
    delta = 0.0;                                        //detuning
    
    ///Time propagation
    tmax = 1.0e0;
    dt = 1.0e-6;        //second
    notimesteps = tmax/dt;
    printf("%d\n", notimesteps);
    
    //store vector
    excited_num = malloc (Nstates*sizeof(unsigned));        //excited number
    y = malloc (Nbasis*sizeof(double));                     //basis
    mismatch = malloc (Natoms * sizeof(double));             //mismatch energies
    atom_x = (double*)malloc(Natoms *sizeof(double*));        //coordinate_x
    atom_y = (double*)malloc(Natoms *sizeof(double*));        //coordintae_y

}

void see_coupling (void){
    
    double rad_ij;
    double Nrad;
    double Vab_l;
    double Vab_s;
    double Gab;
    
    Nrad = 800;                                             //nm
    char filename_coup [80];
    FILE *outcoup;
    
    sprintf(filename_coup, "coupling%datoms_long.txt", Natoms);

    outcoup = fopen(filename_coup, "w");
    
    rad_ij = 0.0;
    while (rad_ij < Nrad){
        
        Gab = rad_decay(rad_ij);
        Vab_l = V_hop (rad_ij);
        Vab_s = V0*exp(-rad_ij/rt);
        rad_ij += 0.1;
        
        fprintf(outcoup, "%.1lf\t%.6e\t%.6e\t%.6e\n", rad_ij, Vab_l, Vab_s, Gab);
    }
    
    fclose(outcoup);
}

void write_data_atom (int run){
    
    sprintf (filename, "%datoms_dim%d_Ga%.2lf_R%.2lf_re%d_inv.dat", Natoms, dim, G0, R, run);
    
    if ((outfile = fopen(filename, "w"))==NULL){
        fprintf (stderr, "error writing file %s", filename);
        exit(-7);
    }
}

void free_alloc_atoms(void){
    
    free(excited_num);
    free(y);
    free(mismatch);
    free(atom_x);
    free(atom_y);
    
}

void init_ground_atom (void){
    
    for (bra = 0; bra < Nstates; bra+=1){
        
        for (ket = 0; ket < Nstates; ket+=1){
            
            y[bra*Nstates + ket] = 0.0;
        }
    }
    
    y[2*(1*Nstates + 1)] = 1.0;
}

void init_inverted_atom (void){
    
    for (bra = 0; bra < Nstates; bra+=1){
        
        for (ket = 0; ket < Nstates; ket+=1){
            
            y[bra*Nstates + ket] = 0.0;
        }
    }
    
    y[Nbasis-2] = 1.0;
}

void init_symmetry_atom(void){
    for (bra = 0; bra < Nstates; bra+=1){
        for (ket = 0; ket < Nstates; ket+=1){
            for (atom_i = 0; atom_i < Natoms; atom_i+=1){
                for (atom_j = 0; atom_j < Natoms; atom_j+=1){
                    if( ((1<<atom_j)==ket) && ((1<<atom_i)==bra))y[2*(bra*Nstates + ket)]=1.0/Natoms;
                }
            }
        }
    }
}

void init_configuration_atoms(void){

    ///initializing global noise
    
    gettimeofday(&tv, NULL);
    nowtime = tv.tv_sec;
    init_genrand(tv.tv_usec);
    
    ///initializing local noise & basis
    for (state=0; state < Nstates; state+=1){
        excited_num[state] = 0;
    }
    write_data_atom(realisation);
    for (atom_i = 0; atom_i<Natoms; atom_i+=1){
        mismatch[atom_i] = ssigma*rand_gauss();
//        fprintf(outfile, "%.6e\n", mismatch[atom_i]);
    }
    
//    init_ground_atom();
    init_inverted_atom();
//    init_symmetry_atom();

///Determining the geometry of atoms arrangement
    twodrandposition(); 
    
}

void twodrandposition (void){
  ///Determining random position of atoms
    atom_i = 0;
    while(atom_i < Natoms){
        atom_x[atom_i]= genrand_real2()*Lx;
        atom_y[atom_i]= genrand_real2()*Ly;
        
        change=0;
        for(atom_j=0; atom_j<atom_i;atom_j++){
            rad=sqrt((atom_x[atom_i]-atom_x[atom_j])*(atom_x[atom_i]-atom_x[atom_j])
                     +(atom_y[atom_i]-atom_y[atom_j])*(atom_y[atom_i]-atom_y[atom_j]));
            
            if(rad < r_min)
                change=1;
        }
        if(change==0)
            atom_i++;
    }

   /*  
    for (atom_i = 0; atom_i < Natoms-1; atom_i+=1){
        rad=sqrt((atom_x[atom_i]-atom_x[atom_i+1])*(atom_x[atom_i]-atom_x[atom_i+1])
                 +(atom_y[atom_i]-atom_y[atom_i+1])*(atom_y[atom_i]-atom_y[atom_i+1]));
        printf("%.2lf, %.2lf \t %.2lf \n", atom_x[atom_i], atom_y[atom_i], rad);
    }
    */

}


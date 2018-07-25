/*
 * =====================================================================================
 *
 *       Filename:  carana_radsup.c
 *
 *    Description:  Kumaha kajadian radiasi super nepika bisa dina dunia nyata
 *
 *        Version:  1.0
 *        Created:  02/23/2018 17:39:10
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wildan Abdussalam (), edun_as@icloud.com
 *   Organization:  Ziza institut
 *
 * =====================================================================================
 */

#include "sirah_radsup.h"

double rand_gauss (void){
    
    static double V1, V2, S;
    static int phase = 0;
    double X;
    
    if(phase == 0) {
        do {
            double U1 = genrand_real2();
            double U2 = genrand_real2();
            
            V1 = 2 * U1 - 1;
            V2 = 2 * U2 - 1;
            S = V1 * V1 + V2 * V2;
        } while(S >= 1 || S == 0);
        
        X = V1 * sqrt(-2 * log(S) / S);
    } else
        
        X = V2 * sqrt(-2 * log(S) / S);
    phase = 1 - phase;
    return X;
}
//status OK
double V_hop (double rij){
    
    double ksi, ksi2, ksi3;
    double total;
    
    ksi=2*M_PI*n_refr*rij/lambda;
	ksi2=ksi*ksi;
	ksi3=ksi*ksi2;
    
    total = 0.5 * 0.75 * G0 * ( cos(ksi)/ksi + sin(ksi)/ksi2 + cos(ksi)/ksi3 );
    total += V0*exp(-rij/rt);
    
	return total;
}

//status OK
double rad_decay(double rij){
    
    double ksi, ksi2, ksi3;
    double total;
    
    ksi=2*M_PI*n_refr*rij/lambda;
	ksi2=ksi*ksi;
	ksi3=ksi*ksi2;
    
    total = 0.5 * 1.50 * ( sin(ksi)/ksi - cos(ksi)/ksi2 + sin(ksi)/ksi3 );
    total *= G0;
    
    if (rij==0)total = G0;
    
	return total;
    
}

//status OK
void counting_excited_state_atoms(void){
    
    for (state = 0; state < Nstates; state+=1){
        for (atom_i = Nphotons; atom_i < Natoms + Nphotons; atom_i+=1){
            if ((1<<atom_i)&state){
                excited_num[state]+=1;
            }
        }
    }
    
}

int
funcatom (double t, const double y[], double f[],
      void *params){

    double sumreal;
    double sumimag;
    double dx, dy, dx2, dy2;
    double Gij;
    double Vij;
    int sigma_up;
    int sigma_down;
    int a,b;

    
    struct ODE_params *p = (struct ODE_params *)params;
    double Delta0 = p->Delta0;
    
    for (bra = 0; bra < Nstates; bra+=1){
        for (ket = 0; ket < Nstates; ket+=1){
            sumimag = 0.0;
            sumreal = 0.0;
///laser detuning            
            sumimag +=  Delta0 * y[2*(bra*Nstates + ket) + 1] * excited_num[ket];
            sumimag -=  Delta0 * y[2*(bra*Nstates + ket) + 1] * excited_num[bra];
            sumreal -= Delta0 * y[2*(bra*Nstates + ket)] * excited_num[ket];
            sumreal += Delta0 * y[2*(bra*Nstates + ket)] * excited_num[bra];
            
            for (atom_i = Nphotons; atom_i < (Nphotons + Natoms); atom_i+=1) {
///local laser detuning
                if (ssigma >0){
                    if (ket & (1<<atom_i)){
                        sumimag += mismatch[atom_i] * y[2*(bra*Nstates + ket)+1];
                        sumreal -= mismatch[atom_i] * y[2*(bra*Nstates + ket)];
                    }

                    if (bra & (1<<atom_i)){
                        sumimag -= mismatch[atom_i] * y[2*(bra*Nstates + ket)+1];
                        sumreal += mismatch[atom_i] * y[2*(bra*Nstates + ket)];
                    }
                }
///Pumping rate                
                if (!((1<<atom_i)&ket)){
                    sumimag -=  R/2.0 * y[2*(bra*Nstates + ket)] ;
                    sumreal -=  R/2.0 * y[2*(bra*Nstates + ket)+1] ;
                }

                if (!((1<<atom_i)&bra)){
                    sumimag -=  R/2.0 * y[2*(bra*Nstates + ket)];
                    sumreal -=  R/2.0 * y[2*(bra*Nstates + ket)+1] ;
                }

                if ((1<<atom_i)&ket && (1<<atom_i)&bra){
                    a = (1<<atom_i)^ket;
                    b = (1<<atom_i)^bra;

                    sumimag += R * y[2*(b * Nstates + a)];
                    sumreal += R * y[2*(b * Nstates + a)+1];
                }

                for (photon = 0; photon < Nphotons; photon+=1){
///atom-cavity coupling
                    if ((1<<photon)&ket){
                        sigma_down = (1<<photon)^ket;
                        if (!((1<<atom_i)&sigma_down)){
                            sigma_up = (1<<atom_i)^sigma_down;
                            sumimag += gcav * y[2*(bra * Nstates + sigma_up)+1];
                            sumreal -= gcav * y[2*(bra * Nstates + sigma_up)];
                        }

                    }

                    if ((1<<atom_i)&ket){
                        sigma_down = (1<<atom_i)^ket;
                        if (!((1<<photon)&sigma_down)){
                            sigma_up = (1<<photon)^sigma_down;
                            sumimag += gcav * y[2*(bra * Nstates + sigma_up)+1];
                            sumreal -= gcav * y[2*(bra * Nstates + sigma_up)];
                        }

                    }

                    if ((1<<photon)&bra){
                        sigma_down = (1<<photon)^bra;
                        if (!((1<<atom_i)&sigma_down)){
                            sigma_up = (1<<atom_i)^sigma_down;
                            sumimag -= gcav * y[2*(sigma_up * Nstates + ket)+1];
                            sumreal += gcav * y[2*(sigma_up * Nstates + ket)];
                        }
                    }

                    if ((1<<atom_i)&bra){
                        sigma_down = (1<<atom_i)^bra;
                        if (!((1<<photon)&sigma_down)){
                            sigma_up = (1<<photon)^sigma_down;
                            sumimag -= gcav * y[2*(sigma_up * Nstates + ket)+1];
                            sumreal += gcav * y[2*(sigma_up * Nstates + ket)];
                        }
                    }

///cavity loss
                    if ((1<<photon)&ket){
                        sumimag -=  kappa/2.0 * y[2*(bra*Nstates + ket)] ;
                        sumreal -=  kappa/2.0 * y[2*(bra*Nstates + ket)+1] ;
                    }

                    if ((1<<photon)&bra){
                        sumimag -=  kappa/2.0 * y[2*(bra*Nstates + ket)];
                        sumreal -=  kappa/2.0 * y[2*(bra*Nstates + ket)+1] ;
                    }

                    if (!((1<<photon)&ket) && !((1<<photon)&bra)){
                        a = (1<<photon)^ket;
                        b = (1<<photon)^bra;

                        sumimag += kappa*y[2*(b * Nstates + a)];
                        sumreal += kappa*y[2*(b * Nstates + a)+1];
                    }
                }

                for (atom_j = Nphotons; atom_j < (Nphotons + Natoms); atom_j+=1){

                  if (dim == 1){
                    rad = abs(atom_i - atom_j);
                  }
                  if (dim == 2){
                    dx = atom_x[atom_j-Nphotons]-atom_x[atom_i-Nphotons];
                    dx2 = dx*dx;
                    dy = atom_y[atom_j-Nphotons]-atom_y[atom_i-Nphotons];
                    dy2 = dy*dy;
                    rad = sqrt(dx2 + dy2);
                  }
///dephasing -- activate this for collective dephasing
                    /*
                    if ( (ket&(1<<atom_i))&&(ket&(1<<atom_j))){
                        sumimag-=gdeph/2.0 *  y[2*(bra*Nstates + ket)];
                        sumreal-=gdeph/2.0 *  y[2*(bra*Nstates + ket)+1];
                    }

                    if ( (bra&(1<<atom_i))&&(bra&(1<<atom_j))){
                        sumimag-=gdeph/2.0 *  y[2*(bra*Nstates + ket)];
                        sumreal-=gdeph/2.0 *  y[2*(bra*Nstates + ket)+1];
                    }

                    if (bra&(1<<atom_i)&&ket&(1<<atom_j)){
                        sumimag+=gdeph * y[2*(bra*Nstates + ket)];
                        sumreal+=gdeph * y[2*(bra*Nstates + ket)+1];
                    }
                    */

                    Gij = rad_decay(rad);

///spontaneous decay
                    if ((1<<atom_j)&ket){
                        sigma_down = (1<<atom_j)^ket;
                        if (!((1<<atom_i)&sigma_down)){
                            sigma_up = (1<<atom_i)^sigma_down;
                            sumimag -= Gij/2.0* y[2*(bra * Nstates + sigma_up)];
                            sumreal -= Gij/2.0* y[2*(bra * Nstates + sigma_up)+1];
                        }
                        
                    }
                      
                    if ((1<<atom_i)&bra){
                        sigma_down = (1<<atom_i)^bra;
                        if (!((1<<atom_j)&sigma_down)){
                            sigma_up = (1<<atom_j)^sigma_down;
                            sumimag -= Gij/2.0 * y[2*(sigma_up * Nstates + ket)];
                            sumreal -= Gij/2.0 * y[2*(sigma_up * Nstates + ket)+1];
                        }
                    }
                      
                    if (!((1<<atom_j)&bra) && !((1<<atom_i)&ket)){
                        a = (1<<atom_j)^bra;
                        b = (1<<atom_i)^ket;
                        sumimag += Gij * y[2*(a * Nstates + b)];
                        sumreal += Gij * y[2*(a * Nstates + b)+1];
                        
                    }
                        
                  if (rad > 0){
                    Vij = V_hop(rad);
                    //printf("%d\t%d\t%.3lf\n", atom_i, atom_j, rad);
                        
///vdW interaction -- activate this for sigma_zz interaction type
                    /*
                    if ( (ket&(1<<atom_i))&&(ket&(1<<atom_j))){
    
                        sumimag += V0/2.0 * y[2*(bra*Nstates + ket)+1]
                        /(rsquare*rsquare*rsquare);
    
                        sumreal -= V0/2.0 * y[2*(bra*Nstates + ket)]
                        /(rsquare*rsquare*rsquare);
                    }
        
                    if ( (bra&(1<<atom_i))&&(bra&(1<<atom_j))){
    
                        sumimag -= V0/2.0 * y[2*(bra*Nstates + ket)+1]
                        /(rsquare*rsquare*rsquare);
    
                        sumreal += V0/2.0 * y[2*(bra*Nstates + ket)]
                        /(rsquare*rsquare*rsquare);
    
                    }
                    */
                        
///forster interactions
                    if ((1<<atom_j)&ket){
                        sigma_down = (1<<atom_j)^ket;
                        if (!((1<<atom_i)&sigma_down)){
                            sigma_up = (1<<atom_i)^sigma_down;
                            
                            sumimag += Vij * y[2*(bra * Nstates + sigma_up)+1];
                            sumreal -= Vij * y[2*(bra * Nstates + sigma_up)];
                            
                        }
                        
                    }
                    
                    if ((1<<atom_i)&bra){
                        sigma_down = (1<<atom_i)^bra;
                        if (!((1<<atom_j)&sigma_down)){
                            sigma_up = (1<<atom_j)^sigma_down;
                            
                            sumimag -= Vij * y[2*(sigma_up * Nstates + ket)+1];
                            sumreal += Vij * y[2*(sigma_up * Nstates + ket)];
                        }
                    }
                  }

                }

            }
        
        f[2*(bra*Nstates + ket)] = sumimag;
        f[2*(bra*Nstates + ket) + 1] = sumreal;
        
        }
    }
    
    return GSL_SUCCESS;
    
}

void evolution_gsl_adaptive_atom (double *detuning){
    
    double expectation;
    double expectation2;
    double *lum;
    double difflum;
    double dt1;
    
    const gsl_odeiv2_step_type * T
    = gsl_odeiv2_step_rk8pd;
    
    write_data_atom(realisation);
    
    gsl_odeiv2_step * s
    = gsl_odeiv2_step_alloc (T, Nbasis);
    gsl_odeiv2_control * c
    = gsl_odeiv2_control_y_new (1e-12, 1e-6);
    gsl_odeiv2_evolve * e
    = gsl_odeiv2_evolve_alloc (Nbasis);
    
    struct ODE_params rod_param={*detuning};//ran26
    
    gsl_odeiv2_system sys = {funcatom, 0, Nbasis, &rod_param};
    
    //    init_basis_density();
    t = 0.0;
    Niter = 1e3;
    lum = malloc(Niter * sizeof(double));
    h = 1e-6;
    dt1 = tmax /Niter;
    //printf ("%d\t%lf\t%lf\n", numiter, tmax, h);
    for (iter=0; iter<Niter ; iter+=1)
    {
        t1 = iter * dt1;
        
        while (t < t1)
        {
            int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);
            
            if (status != GSL_SUCCESS){
                break;
            }
        }
/*  
        printf ("%.3lf\t", t);
        fprintf (outfile, "%.3lf\t", t);

        for (basis = 0; basis < Nbasis; basis+=1){
            printf ("%.6e\t", y[basis]);
            fprintf (outfile, "%.6e\t", y[basis]);
        }
        printf ("\n");
        fprintf (outfile, "\n");
*/
        

        expectation = 0.0;
        expectation2 = 0.0;
        norm = 0.0;
        printf ( "%.3lf\t", t );
        fprintf (outfile, "%.3lf\t", t );
        for (bra = 0; bra < Nstates; bra+=1){

            printf("%.6e\t", y[2*(bra*Nstates + bra)]);
            fprintf(outfile, "%.6e\t", y[2*(bra*Nstates + bra)]);
            expectation += (y[2*(bra*Nstates + bra)])*excited_num[bra];
            norm += (y[2*(bra*Nstates + bra)]);
            lum[iter] =expectation;
            expectation2 += (y[2*(bra*Nstates + bra)])*excited_num[bra]*excited_num[bra];
            
        }
        printf ( "%.6e\n", norm );
        fprintf ( outfile, "\n" );
//        fprintf(outfile, "%.3lf\t%.6e\t%.6e\n", t, expectation, expectation2);
//        printf("%.3lf\t%.6e\t%.6e\n", t, expectation, expectation2);
    }
    
    FILE *outfile2;
    char filename_lum [80];
    
    sprintf(filename_lum, "%dLum_dim%d_Ga%.2lf_R%.2lf_re%d_inv.dat", 
        Natoms, dim, G0, R, realisation);
    
    outfile2 = fopen(filename_lum, "w");
    t = 0;
    for (iter = 1; iter < Niter; iter+=1){
        difflum = fabs(lum[iter -1 ] - lum [iter])/dt1;
        fprintf(outfile2, "%.3lf\t%.6e\n", t, difflum);
        t += dt1;
    }
    
    fclose(outfile2);
    free(lum);
    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
}

void evolution_gsl_nonadaptive_atom (double *detuning){
    
    double expectation;
    double expectation2;
    double y1, y2;
    int m;
    
    const gsl_odeiv2_step_type * T
    = gsl_odeiv2_step_rk8pd;
    
    write_data_atom(realisation);
    
    gsl_odeiv2_step * s
    = gsl_odeiv2_step_alloc (T, Nbasis);
    gsl_odeiv2_control * c
    = gsl_odeiv2_control_y_new (1e-12, 1e-6);
    gsl_odeiv2_evolve * e
    = gsl_odeiv2_evolve_alloc (Nbasis);
    
    struct ODE_params rod_param={*detuning};//ran26
    
    gsl_odeiv2_system sys = {funcatom, 0, Nbasis, &rod_param};
    
    //    init_basis_density();
    t = 0.0;
    h = 1e-3;
    
    expectation = 0.0;
    expectation2 = 0.0;
    for (bra = 0; bra < Nstates; bra+=1){
        expectation += (y[2*(bra*Nstates + bra)])*excited_num[bra];
        expectation2 += (y[2*(bra*Nstates + bra)])*excited_num[bra]*excited_num[bra];
        
    }
    
    
    y1 = expectation;
    
    fprintf(outfile, "%.3lf\t%.6e\t%.6e\n", t, expectation, expectation2);
    printf("%.3lf\t%.6e\t%.6e\n", t, expectation, expectation2);
    m = 1;
    while (t < tmax)
    {
        int status = gsl_odeiv2_evolve_apply_fixed_step (e, c, s, &sys, &t, h, y);
        
        if (status != GSL_SUCCESS){
            break;
        }
        
        //        printf ("%.3lf\t", t);
        //        fprintf (outfile, "%.3lf\t", t);
        //
        //        for (basis = 0; basis < Nbasis; basis+=1){
        //            printf ("%lf\t", y[basis]);
        //            fprintf (outfile, "%lf\t", y[basis]);
        //        }
        //        printf ("\n");
        //        fprintf (outfile, "\n");
        
        expectation = 0.0;
        expectation2 = 0.0;
        for (bra = 0; bra < Nstates; bra+=1){
            expectation += (y[2*(bra*Nstates + bra)])*excited_num[bra];
            expectation2 += (y[2*(bra*Nstates + bra)])*excited_num[bra]*excited_num[bra];
            
        }
        
        if (m%2==0)y1 = expectation;
        else y2 = expectation;
        
        fprintf(outfile, "%.3lf\t%.6e\t%.6e\t%.12e\n", t, expectation, expectation2, fabs(y2-y1)/h);
        printf("%.3lf\t%.6e\t%.6e\t%.12e\n", t, expectation, expectation2, fabs(y2-y1)/h);
        
        m+=1;
    }
    
    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);


}


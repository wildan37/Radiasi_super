/*
 * =====================================================================================
 *
 *       Filename:  inti_radsup.c
 *
 *    Description:  Program inti tina radiasi super Strontium nu kakungkung ku rohang
 *
 *        Version:  1.0
 *        Created:  02/23/2018 17:43:02
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Wildan Abdussalam (), edun_as@icloud.com
 *   Organization:  Ziza Institute
 *
 * =====================================================================================
 */

#include "sirah_radsup.h"

int main(int argc, const char * argv[])
{
    if(argc != 4){
        printf("Incorrect number of args\n");
        exit(-1);
    }
    
    if (sscanf(argv[1],"%d", &Natoms) != 1){
        printf("failed to insert Natoms.\n");
        exit(-2);
    }
    
    if (sscanf(argv[2],"%lf", &sigma) != 1){
        printf("failed to insert mismatch energy.\n");
        exit(-3);
    }
    
    if (sscanf(argv[3],"%d", &realisation) != 1){
        printf("failed to insert realisation.\n");
        exit(-4);
    }
    
//    init_params_atoms();
    init_params_srclock();
    init_configuration_atoms();
//    see_coupling();
    counting_excited_state_atoms();
//    evolution_gsl_nonadaptive_atom(&delta);
    evolution_gsl_adaptive_atom (&delta);
//    checking_rad();
    free_alloc_atoms();
    fclose(outfile);

    return 0;
}

/* Name: LMSE
   Author: Miguel Navascués
   Date: 17 April 2008
   Description: Estimation of demographic parameters from linked microsatellite data
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_sort.h>
#include "main.h"
#include "inputoutput.h"
#include "statistics.h"

int main(int argc, char *argv[]){

    clock_t start;
    int i, j;
    int ret;
    int locus, pair;
    double result_MM[3] = {-99 , 100, -99};
    double result_MCL[3] = {-99 , -99 , -99};
    double result_MMc[3] = {-99 , -99, -99};
    double result_MCLc[3] = {-99 , -99 , -99};
    double result_MCLH[3] = {-99 , -99 , -99};
    double indices[5];
    char *opt_file;
    char log_file[12]="Log.txt";
    char csv[5]=".csv";
    char pop_file[60]="Table_";
    double **table_pairs_pop;
    FILE *foutput;
    FILE *fp;
    FILE *fpop;
    int counter_diff;
    double raggedness;
    int maximum=0;
    //double **results;
    //double **table_stats;
    double correction_factor;
    double correction_difference;
    int numerador,denominador;
    double temp1,temp2;

    //Starts counting time
    start=clock();

    //Opens log file
    flog = fopen( log_file, "wt" );
    if (flog) printf("Opening log file %s \n", log_file);
    else {
        printf("Cannot open log file %s \n", log_file);
        return 1;
    }

    //Writing heading in command promt window
    printf("\n\nL * M * S * E\n\n");
    printf("by Miguel Navascues\n\n");
    printf("Program to obtain maximum composite likelihood estimates of demographic\n");
    printf("parameters theta0, theta1 & tau of a stepwise demographic\n");
    printf("expansion from haplotypic data of linked microsatellite loci.\n\n");
    printf("This program uses the GNU Scientific Library 1.8\n");
    printf("http://www.gnu.org/software/gsl\n\n");
    printf("This program uses maximization algorithm implemented by Mark G. Johnson\n");
    printf("http://www.netlib.org/opt/hooke.c\n\n");
    printf("version of the 17th April 2008\n\n\n\n\n");

    //Writing heading in log file
    fprintf(flog,"L * M * S * E\n\n");
    fprintf(flog,"by Miguel Navascues\n\n");
    fprintf(flog,"Program to obtain maximum likelihood estimates of demographic ");
    fprintf(flog,"parameters theta0, theta1 & tau of a stepwise demographic ");
    fprintf(flog,"expansion from haplotypic data of linked microsatellite loci.\n\n");
    fflush(flog);

    printf("\n\n argc = %d \n\n",argc);

    if (argc>1){
        printf("argv: %s %s\n\n",argv[0],argv[1]);
        opt_file=argv[1];
    }


    //Reads options for analysis
    read_options(opt_file);

    //  First READ of INPUT DATA FILE
    read_input(input_file);

    // if statistics will be calculated, it need the following tables:
    // double results[][] for storing the estimates, ragedness index, -logCL...
    // double table_stats[][] to store stats results
    if (statistics==1){
        results=dmatrix(0,20,0,number_of_pop);
        table_stats=dmatrix(0,12,0,number_of_pop);
    }

    /** Vector, mutation & likelihood matrix
    0-number_of_loci-1:   vector of absolute differences (DELTA)
    number of loci:       sum of absolute differences
    number of loci+1:     estimated number of mutations (might me empty)
    number of loci+2:     sum of squared differences
    number_of_loci+3:     [freq in pop] should be empty
    number_of_loci+4:     control (0=unused row, 1=used row)
    number_of_loci+4+x+1: p(DELTA;j) x=0-299
        where x=(j-table[pair][number_of_loci])/2                            **/
    table_mutations=dmatrix(0,max_num_of_pairs*2,0,number_of_loci+4+300);
    for (i = 0 ; i < max_num_of_pairs*2 ; i++){
        for (j = 0 ; j < number_of_loci+4+300 ; j++){
            table_mutations[i][j]=-9;
        }
    }
    num_of_stored_pairs=0;

    //Calculates relative mutation rates following the chosen approach
    if (option_rates==1){
        rates_opt1(input_file);
        fprintf(flog,"\nRelative mutation rates estimated from allele size variance.\n");
    }
    if (option_rates==2){
        rates_opt2();
        fprintf(flog,"\nMutation rates considered uniform among loci.\n");
    }
    if (option_rates==3) rates_opt3();
    if (option_rates==4){
        rates_opt4(input_file);
        fprintf(flog,"\nRelative mutation rates estimated from homozigosity.\n");
    }


    //prints in log file values of relative mutation rates (multinomial parameters)
    fprintf(flog,"Multinomial parameters are (from relative mutation rates):\n");
	for (locus = 0 ; locus < number_of_loci ; locus++){
		fprintf(flog,"p[locus %d] = %e\n",locus+1,p[locus]);
	}
    fprintf(flog, "\n\n");
    fflush(flog);

    if (option_rates==1 || option_rates==4){
        //prints in log file values of relative mutation rates to locus 1
        fprintf(flog,"Relative mutation rates respect to locus 1 (95 0/0 CI):\n");
        for (locus = 0 ; locus < number_of_loci ; locus++){
            fprintf(flog,"mut[locus %d] = %e(%e-%e)\n",locus+1,mut[locus],mut95[locus][0],mut95[locus][1]);
        }
        fprintf(flog,"total sample size = %d\n",total_sample_size);
        fprintf(flog, "\n\n");
        fflush(flog);
    }

    // OPENS input DATA file again to read it population by population
    fp = fopen(input_file,"r");
    if(fp) printf( "\nOpening input file %s \n", input_file );
    else{
        printf( "Cannot open file %s \n", input_file );
        return 1;
    }

    // READS head of INPUT DATA file
    read_head(fp);

    // OPENS input RESULTS fileto save results on it (right)
    foutput = fopen( output_file, "wt" );
    if( foutput ) printf( "Opening output file %s \n", output_file );
    else{
        printf( "Cannot open output file %s \n", output_file );
        return 1;
    }

    // WRITES head of RESULTS output file, depending on the different method used
    writes_head(foutput);

    /** LOOP for populations STARTS HERE: reads & analysed pop by pop **/
    for (pop = 0 ; pop < number_of_pop ; pop++){

        problem=0;

        //reads input file
        read_1pop(fp);

        printf("\nPop: \"%s\" (%d of %d)\n",pop_name[pop],pop+1,number_of_pop);
        printf("**********************\n\n");

        //sets table to store the "mismatch distribution"
        table_pairs_pop=dmatrix(0,3,0,num_of_pairs[pop]);
        for (i=0;i<3;i++){
            for (j=0;j<num_of_pairs[pop];j++) table_pairs_pop[i][j]=0;
        }

        //sets values of estimates to error flag -9
        for (i=0;i<3;i++){
            result_MM[i] = -9;
            result_MCL[i] = -9;
            result_MCLH[i] = -9;
        }

        //prints population name to results file
        fprintf(foutput,"\"%s\"", pop_name[pop]);

        if (save_mismatch_option==1){
            //sets name of population mismatch distribution output file
            pop_file[0]='\0';
            strcat( pop_file, pop_name[pop] );
            strcat( pop_file, csv );
            //opens population mismatch distribution output file
            fpop = fopen( pop_file, "wt" );
            if( fpop ) printf( "Opening population file %s \n", pop_file );
            else{
                printf( "Cannot open population file %s \n", pop_file );
                return 1;
            }
            //writes head of population mismatch distribution output file
            fprintf(fpop,"\"sum of absolute differences\",\"count\",\"frequency\"\n");
        }

        //calculates the homoplasy factor = # mutations/# differences
        if(method_MM_c==1 || method_MCL_c==1 || method_MCLH==1){
            correction_factor=0.0;
            correction_difference=0.0;
            numerador=0;
            denominador=0;
            maximum=0;
            printf("\nCalculating ML estimates of mutation number\n");
            for (pair = 0 ; pair < num_of_pairs[pop] ; pair++){
                if (delta_vector[pair][number_of_loci+4]==1.0){
                    printf("*");//if ((pair % 50)==0) printf("-");
                    delta_vector[pair][number_of_loci+1] = (double)max_mut_SMM(pair);

                    if (delta_vector[pair][number_of_loci]>maximum)
                        maximum=(int)delta_vector[pair][number_of_loci];
                    correction_difference = correction_difference
                                        + (delta_vector[pair][number_of_loci+1]
                                        * delta_vector[pair][number_of_loci+3])
                                        - (delta_vector[pair][number_of_loci]
                                        * delta_vector[pair][number_of_loci+3]);
                    numerador += (delta_vector[pair][number_of_loci+1]
                                        * delta_vector[pair][number_of_loci+3]);
                    denominador += (delta_vector[pair][number_of_loci]
                                        * delta_vector[pair][number_of_loci+3]);
                }
            }
            correction_factor = (double)numerador / (double)denominador;
            correction_difference = correction_difference / (double)num_of_pairs[pop];
            printf("\n");
            printf("Homoplasy factor = %g\n",correction_factor);
        }else{
            correction_factor=-9;
        }

        //calculates mismatch distributions and saves it to a file

        for (i=0;i<=maximum+1;i++){
            counter_diff=0;
            for (pair = 0 ; pair < num_of_pairs[pop] ; pair++){
                if ((int)delta_vector[pair][number_of_loci]==i) counter_diff++;
            }
            table_pairs_pop[0][i]=counter_diff;
            table_pairs_pop[1][i]=(double)counter_diff/(double)num_of_pairs[pop];
            if (save_mismatch_option==1) fprintf(fpop,"%d,%g,%g\n",i,table_pairs_pop[0][i],table_pairs_pop[1][i]);
        }
        if (save_mismatch_option==1) fflush(fpop);


        //calculates diversity indices
        diversity_index(indices);
        //for (i=0;i<5;i++) printf(" %g ",indices[i]);

        //calculates raggedness index
        raggedness=0;
        for (i=1;i<=maximum+1;i++){
            raggedness=raggedness + pow( table_pairs_pop[1][i] - table_pairs_pop[1][i-1] , 2.0);
        }

        if (statistics==1){
            results[5][pop]=indices[4];
            results[6][pop]=correction_factor;

            results[11][pop]=pop_size[pop];

            results[16][pop]=indices[0];
            results[17][pop]=indices[1];
            results[18][pop]=indices[2];
            results[19][pop]=indices[3];
        }




        //saves to reslts file genetic diversity indices
        fprintf(foutput,",%d",pop_size[pop]);
        fprintf(foutput,",%g",indices[0]);
        fprintf(foutput,",%g",indices[1]);
        fprintf(foutput,",%g",indices[2]);
        fprintf(foutput,",%g",indices[3]);
        fprintf(foutput,",%g",indices[4]);
        fprintf(foutput,",%g",correction_factor);

        //calculates moment estimates of demographic parameters
        printf("\nCalculating moment estimates of demographic parameters\n");
        rogers_1995(result_MM);
        printf("\nEstimates from moments method\n");
        printf("T0=%g\n",result_MM[0]);
        printf("t=%g\n",result_MM[2]);

        fprintf(foutput,",%g", result_MM[0]);
        fprintf(foutput,",%g", result_MM[2]);
        if (method_MM_c==1){
            fprintf(foutput,",%g", result_MM[2] * correction_factor );
            for (i=0;i<2;i++) result_MMc[i]=result_MM[i];
            result_MMc[2]=result_MM[2] * correction_factor;
        }
        fflush(foutput);

        if (statistics==1){
            results[0][pop]=result_MM[0];
            results[1][pop]=result_MM[2];
            results[14][pop]=result_MMc[2];
        }

        option=1;
        if (method_MCL==1){
            printf("\nCalculating ML estimates of demographic parameters (model withOUT homoplasy)\n");
            hooke(3, result_MM, result_MCL, 0.5, 0.0001, 1000);
            if (result_MCL[0]==-9) hooke(3, result_MCL, result_MCL, 0.85, 0.0001, 1000);
            printf("\n");
            printf("\nEstimates from ML method based on model without homoplasy\n");
            printf("T0=%g\n",result_MCL[0]);
            printf("T1=%g\n",result_MCL[1]);
            printf("t=%g\n",result_MCL[2]);

            if (statistics==1){
                results[2][pop]=result_MCL[0];
                results[3][pop]=result_MCL[1];
                results[4][pop]=result_MCL[2];
                if (method_MCL_c==1)
                    results[15][pop]=result_MCL[2]*correction_factor;
            }
            for (i=0;i<3;i++) fprintf(foutput,",%g", result_MCL[i]);
            if (method_MCL_c==1){
                fprintf(foutput,",%g", result_MCL[2]*correction_factor);
                for (i=0;i<2;i++) result_MCLc[i]=result_MCL[i];
                result_MCLc[2]=result_MCL[2] * correction_factor;
            }
            fflush(foutput);
        }


        if (method_MCLH==1){
            option=2;
            printf("\nCalculating MCL estimates of demographic parameters (model WITH homoplasy)\n");
            //if (result_moments_mut[0]==0) result_moments_mut[0]=0.0001;
            //if (result_moments_mut[2]==0) result_moments_mut[2]=0.0001;


            /* NB: empirically (in silico) it seems that using as starting point
            the MCL estimates lowers bias of tau estimate and lowers computation
            time. It also increases a bit the variance of estimates, but mean
            squared error are very similar. So we have decided to use this as
            starting point except for the few cases where MCL fail to converge,
            where moments estimates are used for starting point. */
            if (result_MCL[0]==-9){
                hooke(3, result_MM, result_MCLH, 0.5, 0.0001, 1000);
                if (result_MCLH[0]==-9) hooke(3, result_MM, result_MCLH, 0.85, 0.0001, 1000);
            }else{
                hooke(3, result_MCL, result_MCLH, 0.5, 0.00001, 1000);
                if (result_MCLH[0]==-9) hooke(3, result_MCL, result_MCLH, 0.85, 0.0001, 1000);
            }

            printf("\n");
            printf("\nEstimates from MCL method based on model with homoplasy\n");
            printf("T0=%g\n",result_MCLH[0]);
            printf("T1=%g\n",result_MCLH[1]);
            printf("t=%g\n",result_MCLH[2]);

            for (i=0;i<3;i++) fprintf(foutput,",%g", result_MCLH[i]);
            fflush(foutput);

            if (statistics==1){
                results[7][pop]=result_MCLH[0];
                results[8][pop]=result_MCLH[1];
                results[9][pop]=result_MCLH[2];
            }
        }

        fprintf( foutput , ",%g" , raggedness );
        if (statistics==1){
            results[10][pop]=raggedness;

            option=1;
            results[12][pop]=p_sample_no_H(ref_for_stats);
            option=2;
            results[13][pop]=-9;
            if (method_MCLH==1) results[13][pop]=p_sample(ref_for_stats);
            fprintf( foutput , ",%g" , results[12][pop] );
            fprintf( foutput , ",%g" , results[13][pop] );
        }else{
            if (method_MCL==1){
                option=1;
                temp1=p_sample_no_H(result_MCL);
                fprintf( foutput , ",%g" , temp1 );
            }
            if (method_MCLH==1){
                option=2;
                temp2=p_sample(result_MCLH);
                fprintf( foutput , ",%g" , temp2 );
            }
        }

        printf("\nPop: \"%s\" (%d of %d) finished\n",pop_name[pop],pop+1,number_of_pop);
        printf("*******************************\n\n");

        free_imatrix(data,0,pop_size[pop],0,number_of_loci+1);
        free_dmatrix(delta_vector,0,num_of_pairs[pop],0,number_of_loci+5);
        //free_dmatrix(table_pairs_pop,0,2,0,?????);

        if (problem==0) fprintf(flog, "Analyses on population \"%s\" were done without detected problems\n", pop_name[pop]);
        else fprintf(flog, "Problems arised on the analysis of population population \"%s\"\n", pop_name[pop]);
        fflush(flog);

        fprintf(foutput,"\n");

        if (save_mismatch_option==1){
            ret = fclose(fpop);
            if( ret==0 ) printf( "File %s closed\n", pop_file );
            else printf( "\nCannot close file %s \n", pop_file );
        }

    }
/** LOOP for populations ENDS HERE: reads & analysed pop by pop **/

    if (statistics==1) sim_stats(foutput);

    printf("Computed in %f s\n\n\n\n\n\n\n", (clock()-start)/(double)CLOCKS_PER_SEC );
    fprintf(flog,"\n\n\nComputed in %f s\n", (clock()-start)/(double)CLOCKS_PER_SEC );

    //system("PAUSE");

    return 0;
}

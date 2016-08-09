#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "inputoutput.h"

int read_options(char *opt_file){

    FILE *fp;
    //char opt_file[50]="options.txt";
    int number_of_options,option,i;

    if (opt_file[0]=='\0'){
        printf("\nEnter the name of the option file\n");
        scanf("%s",opt_file);
        fprintf(flog, "\nOption file: %s \n", opt_file );
    }

    //Opens input file with options
    fp = fopen( opt_file, "r" );
    if( fp ) printf( "Opening option file %s \n", opt_file );
    else{
        printf( "Cannot open option file %s \n", opt_file );
        return 1;
    }

    //Default options
    method_MM_c=0;
    method_MCL=0;
    method_MCL_c=0;
    method_MCLH=1;
    statistics=0;
    ref_for_stats[0]=-9;
    ref_for_stats[1]=-9;
    ref_for_stats[2]=-9;
    ref_for_stats[3]=-9;
    ref_for_stats[5]=-9;
    ref_for_stats[6]=-9;
    //input_file="data.txt";
    //output_file="Results.csv";
    option_rates=1;
    //input_rates="rates.txt";
    q=0.5;
    save_mismatch_option=0;




    //Reading options
    option=0;
    number_of_options = 20;
    bootstrap=0;
    for ( i=1 ; i<=number_of_options ; i++ ){
        fscanf(fp,"OPTION %d:",&option);
        if (option==1)  fscanf( fp, "%d \n", &method_MM_c);
        if (option==2)  fscanf( fp, "%d \n", &method_MCL);
        if (option==21) fscanf( fp, "%d \n", &method_MCL_c);
        if (option==3)  fscanf( fp, "%d \n", &method_MCLH);
        if (option==4)  fscanf( fp, "%d \n", &statistics);
        if (option==41) fscanf( fp, "%lf\n", &ref_for_stats[0]);
        if (option==42) fscanf( fp, "%lf\n", &ref_for_stats[1]);
        if (option==43) fscanf( fp, "%lf\n", &ref_for_stats[2]);
        if (option==44) fscanf( fp, "%lf\n", &ref_for_stats[3]);
        if (option==45) fscanf( fp, "%lf\n", &ref_for_stats[5]);
        if (option==46) fscanf( fp, "%lf\n", &ref_for_stats[6]);
        if (option==5)  fscanf( fp, "%s \n", input_file);
        if (option==6)  fscanf( fp, "%s \n", output_file);
        if (option==7)  fscanf( fp, "%d \n", &option_rates);
        if (option==71) fscanf( fp, "%s \n", input_rates);
        if (option==8)  fscanf( fp, "%lf\n", &q);
        if (option==9)  fscanf( fp, "%d \n", &save_mismatch_option);
        if (option==10) fscanf( fp, "%d \n", &bootstrap);
    }

    if (method_MM_c!=1 && method_MM_c!=0) method_MM_c=0;
    if (method_MCL!=1 && method_MCL!=0) method_MCL=0;
    if (method_MCL_c!=1 && method_MCL_c!=0) method_MCL_c=0;
    if (method_MCL_c==1) method_MCL=1;
    if (method_MCLH!=1 && method_MCLH!=0) method_MCLH=1;
    if (method_MCLH==1) method_MCL=1;
    if (statistics!=1 && statistics!=0) statistics=0;
    if (option_rates!=1 && option_rates!=2 && option_rates!=3 && option_rates!=4) option_rates=1;
    if (save_mismatch_option!=1 && save_mismatch_option!=0) save_mismatch_option=0;
    if ( q<=0 || q>=1 ) q=0.5;

    return 0;

}

void read_input(char input_file[50]){

    int pop, ind, locus;
    int ret;
    FILE* fp;
    int trash;

    fp = fopen(input_file,"r");
    if( fp ) printf( "\nOpening input file %s \n", input_file );
    else{
        printf( "Cannot open file %s \n", input_file );
    }

    fprintf(flog,"\nInput file contains:\n");

    fscanf( fp, "%d %d\n", &number_of_pop, &number_of_loci );

    if (number_of_pop==1) printf("\nFile contains 1 population\n");
    else printf("\nFile contains %d populations\n",number_of_pop);

    if (number_of_loci==1) printf("\nIndividuals are typed for 1 locus\n\n");
    else printf("\nIndividuals are typed for %d loci\n\n",number_of_loci);

    fprintf(flog,"%d populations\nindividuals typed for %d loci\n",number_of_pop,number_of_loci);
    fflush(flog);

    pop_size=ivector(0,number_of_pop);
    num_of_pairs=ivector(0,number_of_pop);
    pop_name=cmatrix(0,number_of_pop,0,50);
    p=dvector(0,number_of_loci);
    mut=dvector(0,number_of_loci);
    mut95=dmatrix(0,number_of_loci,0,2);

    total_sample_size=0;
    total_num_of_pairs=0;
    max_num_of_pairs=0;
    for (pop=0 ; pop<number_of_pop ; pop++){
        //printf("pop %d\n",pop);
        fscanf( fp, "%d %s \n", &pop_size[pop], pop_name[pop]);

        num_of_pairs[pop] = pop_size[pop] * ( pop_size[pop]-1 ) / 2;

        printf("Population %d is called \"%s\" and has sample size of %d individuals\n", pop+1, pop_name[pop], pop_size[pop]);
        printf("(this makes %d pairwise haplotypes comparisons)\n\n", num_of_pairs[pop]);

        fprintf(flog, "Population \"%s\" has %d individuals (%d pairwise comparisons)\n", pop_name[pop],pop_size[pop],num_of_pairs[pop]);
        fflush(flog);

        total_sample_size += pop_size[pop];
        total_num_of_pairs += num_of_pairs[pop];

        if(num_of_pairs[pop]>max_num_of_pairs){
            max_num_of_pairs=num_of_pairs[pop];
        }

        for (ind=0 ; ind<pop_size[pop] ; ind++){
            for (locus=0 ; locus<number_of_loci ; locus++) fscanf(fp,"%d",&trash);
            fscanf(fp,"\n");
        }
    }
    fprintf(flog, "\n");

    printf("max_num_of_pairs=%d\n",max_num_of_pairs);
    if(max_num_of_pairs>100000){
        max_num_of_pairs=100000;
        printf("max_num_of_pairs=%d\n",max_num_of_pairs);
    }




    printf("\nTotal sample size is %d individuals\n\n", total_sample_size);

    ret = fclose(fp);

    if( ret==0 )
        printf( "File %s closed\n", input_file );
    else{
        printf( "\nCannot close file %s \n", input_file );
    }
}

void read_head(FILE *fp){

    fscanf( fp, "%d %d\n", &number_of_pop, &number_of_loci );
    //printf("reading head\n");

}

void writes_head(FILE *foutput){

    fprintf(foutput,"\"population\"");

    fprintf(foutput,",\"sample size\"");
    fprintf(foutput,",\"n\"");
    fprintf(foutput,",\"ne\"");
    fprintf(foutput,",\"He\"");
    fprintf(foutput,",\"D2sh\"");
    fprintf(foutput,",\"pi\"");
    fprintf(foutput,",\"homoplasy factor\"");

    fprintf(foutput,",\"Moments theta0\"");
    fprintf(foutput,",\"Moments tau\"");
    if (method_MM_c==1) fprintf(foutput,",\"Moments(c) tau\"");
    if (method_MCL==1){
        fprintf(foutput,",\"C Likelihood theta0\"");
        fprintf(foutput,",\"C Likelihood theta1\"");
        fprintf(foutput,",\"C Likelihood tau\"");
        if (method_MCL_c==1)
            fprintf(foutput,",\"C Likelihood(c) tau\"");
    }
    if (method_MCLH==1){
        fprintf(foutput,",\"C Likelihood(h) theta0\"");
        fprintf(foutput,",\"C Likelihood(h) theta1\"");
        fprintf(foutput,",\"C Likelihood(h) tau\"");
    }
    fprintf(foutput,",\"RI\"");
    fprintf(foutput,",\"-log[CL]\"");
    fprintf(foutput,",\"-log[CL(h)]\"");
    fprintf(foutput,"\n");


}


void read_1pop(FILE *fp){

    int ind,locus,pair;
    int pair_included, mismatch;
    int i, j, k;
    double sum=0.0;
    double SSD=0.0;
    double *temp_vector;

    fscanf( fp, "%d %s \n", &pop_size[pop], pop_name[pop]);

    /** Matrix that contains the raw data **/
    data=imatrix(0,pop_size[pop],0,number_of_loci+1);
    for (i=0 ; i<pop_size[pop] ; i++){
        for (j=0 ; j<number_of_loci+1 ; j++){
            data[i][j]=-9;
        }
    }

    /** Matrix that contains the haplotype pairwise comparisons
    0-number_of_loci-1: vector of differences
    number_of_loci:     sum of absolute differences
    number_of_loci+1:   estimated number of mutations
    number_of_loci+2:   sum of squared differences
    number_of_loci+3:   frequency (count) in pop
    number_of_loci+4:   control (0=unused row, 1=used row) **/

    printf("num_of_pairs[pop]=%d\n",num_of_pairs[pop]);
    if (num_of_pairs[pop]>100000){
        num_of_pairs[pop]=100000;
        printf("num_of_pairs[pop]=%d\n",num_of_pairs[pop]);
    }

    delta_vector=dmatrix(0,num_of_pairs[pop],0,number_of_loci+5);
    for (i=0 ; i<num_of_pairs[pop] ; i++){
        for (j=0 ; j<number_of_loci+5 ; j++){
            delta_vector[i][j]=-9;
        }
    }

    temp_vector=dvector(0,number_of_loci+3);

    for (ind=0 ; ind<pop_size[pop] ; ind++){
        for (locus=0 ; locus<number_of_loci ; locus++){
            fscanf(fp,"%d",&data[ind][locus]);
        }
        fscanf(fp,"\n");
    }

    pair_included=0;
    for (i = 0 ; i < (pop_size[pop] - 1) ; i++){
        for (j = i+1 ; j < pop_size[pop] ; j++){

            //loop to set the values of the vector describing the pairwise comparison
            sum = 0.0;
            SSD = 0.0;
            for (locus=0 ; locus<number_of_loci ; locus++){
                temp_vector[locus]=(double)abs(data[i][locus]-data[j][locus]);
                sum += temp_vector[locus];
                SSD += temp_vector[locus]*temp_vector[locus];
            }
            temp_vector[number_of_loci] = sum;
            temp_vector[number_of_loci+2] = SSD;

            //loop to include the vector in the table if it is new or to increase its freq it is is already present
            for ( k=0 ; k<num_of_pairs[pop] ; k++){
                if (delta_vector[k][number_of_loci+4]!=-9){
                    mismatch=0;
                    if (temp_vector[number_of_loci]==delta_vector[k][number_of_loci]){
                        if (temp_vector[number_of_loci+2]==delta_vector[k][number_of_loci+2]){
                            for (locus=0;locus<number_of_loci;locus++)
                                if (temp_vector[locus]!=delta_vector[k][locus]) mismatch=1;
                        }else mismatch=1;
                    }else mismatch=1;


                    if (!mismatch){
                        delta_vector[k][number_of_loci+3]++;
                        //pair_included++;
                        k=num_of_pairs[pop];
                    }
                }else{
                    for (locus=0 ; locus<number_of_loci+1 ; locus++){
                        delta_vector[k][locus]=temp_vector[locus];
                    }
                    delta_vector[k][number_of_loci+2]=temp_vector[number_of_loci+2];
                    delta_vector[k][number_of_loci+3]=1.0;
                    delta_vector[k][number_of_loci+4]=1.0;
                    pair_included++;
                    k=num_of_pairs[pop];
                }
            }
        }
    }
    num_of_pairs[pop]=pair_included;

    /*printf("\n\n");
    for ( k=0 ; k<num_of_pairs[pop] ; k++){
        for (locus=0 ; locus<number_of_loci+5 ; locus++){
            printf(" %g ", delta_vector[k][locus]);
        }
        printf("\n");
    }
    printf("\n"); */

    //add new vectors to table_muations
    if (num_of_stored_pairs<max_num_of_pairs*2){
        for ( k=0 ; k<num_of_pairs[pop] ; k++){
            if (delta_vector[k][number_of_loci+4]==1.0){
                for (pair=0 ; pair<=num_of_stored_pairs ; pair++){
                    if (pair==num_of_stored_pairs && table_mutations[pair][number_of_loci+4]==-9){
                        for (locus=0 ; locus<number_of_loci+1 ; locus++)
                            table_mutations[pair][locus]=delta_vector[k][locus];
                        table_mutations[pair][number_of_loci+2]=delta_vector[k][number_of_loci+2];
                        table_mutations[pair][number_of_loci+4]=1.0;
                        num_of_stored_pairs++;
                        pair=num_of_stored_pairs+1;
                        //printf("STORED PAIRS: %d\n",num_of_stored_pairs);
                        if (num_of_stored_pairs==max_num_of_pairs*2){
                            k=num_of_pairs[pop];
                        }
                    }else{
                        mismatch=0;
                        if (table_mutations[pair][number_of_loci]==delta_vector[k][number_of_loci]){
                            if (table_mutations[pair][number_of_loci+2]==delta_vector[k][number_of_loci+2]){
                                for (locus=0;locus<number_of_loci;locus++)
                                    if (table_mutations[pair][locus]!=delta_vector[k][locus]) mismatch=1;
                            }else mismatch=1;
                        }else mismatch=1;

                        if (mismatch==0) pair=num_of_stored_pairs+1;
                    }
                }
            }
        }
    }

    /*printf("\n\n");
    for ( pair=0 ; pair<num_of_stored_pairs ; pair++){
        for (locus=0 ; locus<number_of_loci+30 ; locus++){
            printf(" %g ", table_mutations[pair][locus]);
        }
        printf("\n\n\n");
    }
    printf("\n");*/

    free_dvector(temp_vector,0,number_of_loci+3);


}

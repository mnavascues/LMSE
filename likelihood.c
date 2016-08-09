#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <math.h>
//#include <string.h>
//#include <time.h>
#include <assert.h>
#include "likelihood.h"
//#include "memo.h"

#define MAXSIZE 100
#define MAXLOCI 500

/****  function for calculating the factorial of a number ****/
/****  returns the value n! as a floating-point number    ****/
double factorial(int M){

   int i;
   double value=1;
   static float a[12]={1.0,1.0,2.0,6.0,24.0,120.0,720.0,5040.0,40320.0,362880.0,3628800.0,39916800.0};

	if (M < 0){
		printf( "Error: factorial of a negative number\n" );
		return 0;
	}

   if (M<12) value=a[M];
   else{
		value=a[11];
	    for(i=12; i<=M; i++){
			value = value*i;
			// printf("factorial %i %g\n", i, value);
		}
   }
   //printf("factorial of %i = %g\n", M, value);
   //system("PAUSE");
   return value;
}





/**** calculates the probability for getting a number of repeats differences (delta)
given a number of mutations (k) for a microsatellite loci evolving in a strict
stepwise mutation model directionally biased, being 'q' the probability for increasing
number of repeats and '1-q' the probability for decreasing the number of repeats ****/
double p_delta_k(int delta, int k){
    double value;
    double binom1, binom2;

    if (delta == 0){
        if ( (k%2) == 0){
            value = gsl_ran_binomial_pdf(k/2,q,k);
        }else{
            value = 0;
        }
    }else if (delta>0){
        if ( ( (k+delta) % 2 ) == 0){
            binom1 = gsl_ran_binomial_pdf( (k+delta)/2 , q , k);
            binom2 = gsl_ran_binomial_pdf( (k-delta)/2 , q , k);
            value= binom1 + binom2;
        }else{
            value = 0;
        }
    }else{
        printf( "Error: negative delta\n" );
        value = -1;
    }

    if (value < 0){
        printf( "Error: negative probability P(delta;k)\n" );
        printf("P(delta;k) : P(%i;%i) = %g\n", delta,k,value);
        system("PAUSE");
    }

    return value;
}


/**** calculates the probability of getting the vector DELTA that contains the number of
repeat differences on every loci (delta_1, delta_2, ..., delta_L) given a vector K that
contains the number of mutation occurring on every loci (k_1, k_2, ..., k_L) ****/
double p_DELTA_K(int DELTA[], int K[]){
    int i;
    double value = 1;
    for (i=0; i < number_of_loci ; i++){
        value = value * p_delta_k(DELTA[i], K[i]);
    }
    if (value < 0) printf( "Error: negative probability P(DELTA;K)\n" );

    //printf("P(DELTA;K) : P([");
    //for (i=0; i < number_of_loci ; i++)
    //{
    //    printf("%i ",DELTA[i]);
    //}
    //printf("];[",DELTA[i]);
    //for (i=0; i < number_of_loci ; i++)
    //{
    //    printf("%i ",K[i]);
    //}
    //printf("]) = %g\n",value);
    //system("PAUSE");
    return value;
}


/** this function uses recursivity to go through all possible combinations ok j
mutations in distributed on the vector K **/
void combination(int DELTA[], int maxK[], int K[], int counter, int j, double *pointer_sum)
{
    int i, ii, sum, suma;
    double probability_of_combination;
    //double value_p_DELTA_K;
    //double control;

    sum=0;
    //printf("K:");
    for (i=0; i<number_of_loci; i++){
        sum=sum+K[i];
        //printf(" %i",K[i]);
    }
    //printf("; SUM: %i; counter: %i\n", sum, counter);
    if (counter < number_of_loci-1){
       for (i=DELTA[counter]; i<=maxK[counter]; i=i+2){
           K[counter]=i;
           suma=0;
           for (ii=0; ii<number_of_loci; ii++) suma=suma+K[ii];

           // control=p_DELTA_K(DELTA, K) * p_K(K);
           // if ( (suma>j) || (*pointer_maximum_probability_of_combination>control*10) )
           if (suma>j){
              // printf("sum vector K(%i) > j(%i)!!!\n",suma,j);
              goto break_loop;
           }
           /** this is to break the loop when sum of vector K has reached higher
           jalue than j and not to loose time in combination that will never
           give combination compatible with j**/

           combination(DELTA, maxK, K, counter+1, j, pointer_sum);
       }
       break_loop: K[counter]=DELTA[counter];
       /* this resets vector K to initial values */
    }else{
        if (sum<=j){
           K[counter]=DELTA[counter]+j-sum;
           //printf("calculating probability of vector K:");
           //for (i=0; i<number_of_loci; i++)
           //{
           //    printf("%i ",K[i]);
           //}
           //printf("\n");
           //printf("                            y DELTA:");
           //for (i=0; i<number_of_loci; i++)
           //{
           //    printf("%i ",DELTA[i]);
           //}
           //printf("\n");
           suma=0;
           for (i=0; i<number_of_loci; i++) suma=suma+K[i];
           if (suma!=j) printf("ERROR!!! sum(%i)<>j(%i)\n",suma,j);

           probability_of_combination = p_DELTA_K(DELTA, K) * gsl_ran_multinomial_pdf(number_of_loci,p,K);

           *pointer_sum = *pointer_sum + probability_of_combination;
           K[counter]=DELTA[counter];
        }
    }
}


/**** calculates the probability of getting vector DELTA that contains the
number of repet differences on every loci (delta_1, delta_", ..., delta_L) given
the number of mutations, j. For that it integrates over all possible vectors K
that contains how the j mutations are distributed among the L loci ****/
double p_DELTA_j(int DELTA[], int j)
{
    int l, ll;
    int K[MAXLOCI];
    int maxK[MAXLOCI];
    double value = 0;
    double *pointer_sum;
	int pair, locus;
	int pareja;
	int mismatch;
	int x;
	double sum,SSD;

    sum=0;
    SSD=0;
    for (locus=0 ; locus<number_of_loci ; locus++){
        sum += DELTA[locus];
        SSD += DELTA[locus]*DELTA[locus];
    }

    /*printf("DELTA\n");
    for (locus=0 ; locus<number_of_loci+3 ; locus++){
        printf(" %d ",temp_vector[locus]);
    }
    printf("\n");*/

    x=(j-sum)/2;

    //printf("\nX=%d\n",x);
    //printf("J=%d\n",j);
    //printf("sum=%d\n",sum);

	if (x>=299) printf("X%d",x);
	if (x<299){
		for (pair = 0 ; pair < num_of_stored_pairs ; pair++){
            mismatch=0;
            if (sum==table_mutations[pair][number_of_loci]){
                if (SSD==table_mutations[pair][number_of_loci+2]){
                    for (locus=0;locus<number_of_loci;locus++)
                        if (DELTA[locus]!=table_mutations[pair][locus]) mismatch=1;
                }else mismatch=1;
            }else mismatch=1;

            if (mismatch==0){
                pareja = pair;
                pair = num_of_stored_pairs+1;
                //printf("match!!!");
            }
		}
	}

    if (table_mutations[pareja][number_of_loci+4+x+1]!=-9){
        value = table_mutations[pareja][number_of_loci+4+x+1];
        //printf("|");
        return value;
    }

    // printf("INITIAL VALUE = %g\n", value);
    pointer_sum = &value;

    for (l=0; l<number_of_loci; l++) // This calculates the maximum number of mutation
    {                                // at each locus, given the minumum number of mutation
       maxK[l]=j;                    // present in the other loci
       for (ll=0; ll<number_of_loci; ll++){
           if (l!=ll) maxK[l]=maxK[l]-DELTA[ll];
       }
    // printf("max # mut on l %i=%i\n", l, maxK[l]);
    }

    for (l=0; l<number_of_loci; l++){
       K[l]=DELTA[l];
    }

    combination(DELTA, maxK, K,  0, j, pointer_sum);


    if (value < 0){
        printf( "Error: negative probability P(DELTA;j)\n" );
    }
    //printf("P(DELTA;j) : P([");
    //for (l=0; l<number_of_loci; l++)
    //{
    //   printf("%i ", DELTA[l]);
    //}
    //printf("];%i) = %g\n", j,value);
    //system("PAUSE");

   	if (j < 299 && mismatch==0){
        table_mutations[pareja][number_of_loci+4+x+1]=value;
        printf("+");
   	}

	//printf("Calculated value = %g\n",value);

	//printf("delta: ");
	//for (locus = 0 ; locus < number_of_loci ; locus++)
	//{
	//	printf(" %d",DELTA[locus]);
	//}
	//printf("    j: %d\n",j);

	//system("PAUSE");

    return value;
}


/** computes the probability that j mutations occurred in an scenario of constant
populations size with parameter theta = 2*N*mu, where N is the effective pop size
and mu is the mutation rate **/
double p_j_theta(int j, double theta){
       double value;

       value = j * log(theta) - (j+1) * log(theta+1);

       value = exp(value);

       //value = pow(theta,j) / pow(theta+1,j+1);

       //printf("P(j;theta)=%g\n",value);
       if (value < 0) printf( "Error: negative probability P(j;theta)\n" );
       return value;
}


/** computes the probability that j mutations occurred in an scenario of sudden
population growth with parameters theta1 = 2*N1*mu, theta0 = 2*N0*mu, tau=2*t*mu,
where N0 is the effective pop size before the expansion, N1 after the expansion,
mu is the mutation rate and t is the time of the stepwise expansion **/
double p_j_theta_theta_tau(int j, double theta0, double theta1, double tau)
{
    double value=0;
    double sum;
    double sumando;
    double A, B, C, D;
    int i;

    sum=0;
    for (i=0; i<=j; i++){
        assert(j>=0);
        assert(theta0>0);
        assert(theta1>0);
        assert(tau>0);

        //printf("pow(%g,%d) :::::::::::::: %g\n", tau, i, pow(tau,i) );

        if (i==0){
            D = 1;
        }else{
            D = pow(tau,i);
        }
        C = factorial(i);
        //C = factrl(i);
        A = D / C;
        B = p_j_theta(j-i,theta0) - p_j_theta(j-i,theta1);
        sumando = A * B;
        sum = sum + A*B;

           /*if ( (sumando < 0) && (i == 0) ) {
			   printf( "sum = %g\n" , sum );
			   printf( "sumando = %g\n" , sumando );
			   printf( "tau = %g\n" , tau );
			   printf( "D = %f\n" , D );
			   printf( "C = %f\n" , C );
			   printf( "A = %f\n" , A );
			   printf( "B = %f\n" , B );
			   printf( "pow(tau,i) = %g \n" , pow(tau,i) );
			   printf( "factorial(i) = %g \n" , gsl_sf_fact(i) );
			   printf( "p_j_theta(j-i,theta0) = %g\n" , p_j_theta(j-i,theta0) );
			   printf( "p_j_theta(j-i,theta1) = %g\n" , p_j_theta(j-i,theta1) );
			   printf( "i = %d j = %d\n", i, j);
			   system("pause");
           }*/



    }

    value = p_j_theta(j,theta1) + exp( -tau * ((theta1+1)/theta1) ) * sum;

    //if (sum < 0) printf ("SUM<0\n");

       //printf("P(j;theta0,theta1,tau)=%g\n",value);
    if (value < 0){
        //printf("p(j;theta) : p(%d;%g)=%g\n",j,theta1, p_j_theta(j,theta1) );
        //printf("exp=%g\n", exp( -tau * ((theta1+1)/theta1) ) );
        //printf("sum=%g\n", sum );
        printf("WARNING: negative probability, look at the log file\n");
        fprintf(flog,"P(j;T0,T1,t) : P(%d;%g,%g,%g) = %g\n",j,theta0,theta1,tau,value);
        fprintf(flog,"independent checking of the value of this probability is strongly recommended\n");
        fprintf(flog,"L*M*S*E has considered this probability equal to ZERO\n");
        fflush(flog);
        problem = 1;
        //printf("P(j;T0,T1,t) : P(%d;%g,%g,%g) = 0\n",j,theta0,theta1,tau);
        //system("PAUSE");
    }
       //printf("P(j;T0,T1,t) : P(%i;%g,%g,%g) = %g\n",j,theta0,theta1,tau,value);

    /**OJO**/
    if (value < 0) value=0.0;

    assert(value >=0 );
    return value;
}










double p_DELTA_theta_theta_tau(int DELTA[], double theta0, double theta1, double tau)
{
       int l, j;
       int sum;
       double increase = 0;
       double value = 0;

       sum=0;
       for (l=0; l<number_of_loci; l++) sum=sum + DELTA[l];

       j=sum;

       // use 10000000 for higher precission but more time consuming
	   while ( value <= increase*100000){
			//printf(".................................................j=%i\n",j);
			//printf("P(j;theta)=%g\n",p_j_theta_theta_tau(j,theta0,theta1,tau));
			//printf("P(DELTA;j)=%g\n",p_DELTA_j(DELTA,number_of_loci,j,p,q));

			increase = p_j_theta_theta_tau(j,theta0,theta1,tau) * p_DELTA_j(DELTA,j);

			//printf("increase=%g\n",increase);

			value = value + increase;

			//printf("value=%g\n",value);
			//printf("j=%d\n",j);

			j=j+2;
       }

       if (value < 0) printf( "Error: negative probability P(DELTA;theta0,theta1,tau)\n" );

       //printf("P(DELTA;T0,T1,t) : P([");
       //for (l=0; l<number_of_loci; l++)
       //{
       //    printf("%i ", DELTA[l]);
       //}
       //printf("];%g,%g,%g) = %g\n",theta0,theta1,tau,value);
       //system("PAUSE");

       return value;
}










double p_sample(double *punto)
{
	double theta0, theta1, tau;
	int DELTA[MAXLOCI];
	int pair,locus;
	double value = 0;

    //printf("L");

	theta0=punto[0];
	theta1=punto[1];
	tau=punto[2];

	if ( theta0 <= 0 || theta1 <= 0 || tau <= 0){
		return 99999999;
	}else{
       for ( pair = 0 ; pair < num_of_pairs[pop] ; pair++){
           // printf("\nPAIR: %i\n",pair);
           if (delta_vector[pair][number_of_loci+4]){
               // printf("DELTA:");
               for (locus = 0 ; locus < number_of_loci ; locus++){
                   DELTA[locus]=(int)delta_vector[pair][locus];
                   //printf("%i  ",DELTA[locus]);
               }
               //printf("\n");
               value = value + log(p_DELTA_theta_theta_tau(DELTA,theta0,theta1,tau)) * delta_vector[pair][number_of_loci+3];
           }
       }
       value=value*(-1);

	   //printf("\n");

	   return value;
	}
}



double p_sample_no_H(double *punto){

	double theta0, theta1, tau;
	int pair;
	//double A,B;
	int j;
	double value = 0.0;
    double intermediario;

	theta0 = punto[0];
	theta1 = punto[1];
	tau = punto[2];

	if ( (theta0 <= 0) || (theta1 <= 0) || (tau <= 0) ){
		return 99999999;
	}else{
       for ( pair = 0 ; pair < num_of_pairs[pop] ; pair++){
            //printf("\nPAIR: %i\n",pair);
           if (delta_vector[pair][number_of_loci+4]!=1){
               pair=num_of_pairs[pop];
           }else if (delta_vector[pair][number_of_loci+4]==1){
				j = (int)delta_vector[pair][number_of_loci];
				intermediario=p_j_theta_theta_tau(j,theta0,theta1,tau);
				if (intermediario==0){
				    printf("\n\nSome error in p(j,theta0,theta1,tau)\n\n");
				    return 99999999;
                }
                value = value + log(intermediario) * delta_vector[pair][number_of_loci+3];

				//printf("j=%d log-lik=%g \n", j, log(p_j_theta_theta_tau(j,theta0,theta1,tau)));
				//printf("* %d\n",delta_vector[population][pair][number_of_loci+1]);
 				//printf("A=%g B=%g \n", A, B);
				//printf("value=%g A*B=%g \n", value, A*B);

          }
       }
       value=value*(-1.0);

	   //printf("\n");

	   return value;
	}
}


int max_mut_SMM(int pair){

    int j,pareja,par;
    int locus;
    int *DELTA;
    double oldprob, newprob;
    int end=0;
    int match;

    DELTA=ivector(0,number_of_loci);

    //printf("search4max");

    for (locus = 0 ; locus < number_of_loci ; locus++){
        DELTA[locus]=(int)delta_vector[pair][locus];
    }

    for (pareja = 0 ; pareja < num_of_stored_pairs ; pareja++){
        match=1;
        for (locus = 0 ; locus < number_of_loci ; locus++){
            if (DELTA[locus]!=table_mutations[pareja][locus]) match=0;
        }
        if (match==1){
            if (table_mutations[pareja][number_of_loci+1]>=0){
                //printf("useful\n");
                return table_mutations[pareja][number_of_loci+1];
            }else{
                par=pareja;
                pareja=num_of_stored_pairs;
            }
        }
    }

    oldprob = 0.0;
    j=(int)delta_vector[pair][number_of_loci];
    while (end==0){
        //printf("J%d\n",j);
        newprob = p_DELTA_j(DELTA,j);
        if (newprob < oldprob){
            end=1;
            table_mutations[par][number_of_loci+1]=j-2;

            //printf("maxfound: J=%d \n",j-2);
            //system("PAUSE");
            free_ivector(DELTA,0,number_of_loci);
            return j-2;
        }
        oldprob = newprob;
        j=j+2;
    }
    free_ivector(DELTA,0,number_of_loci);

    return -9;
}



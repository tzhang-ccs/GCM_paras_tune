#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define LOGFILE

int NDIM;                    //number of dimensions
int NUM_SIM; 
int UPDATE;                  //update sigal
int NMAX;                    //maximum of iterations
int iterator;                //iterator number
double tolhist;              //tol for stop
double pw;                   //penalty weight
FILE *FP_Log;                //log file pointer
//char benchmark_path[200];    //benchmark path
char benchmark_case[200];     //benchmark executable file
char benchmark_res[200];      //benchmark result file
char constraint_case[200];     //constraint executable file
char constraint_res1[200];      //constraint result file
char constraint_res2[200];      //constraint result file
double **RANGE_Paras;        //paramers range matrix
double **INIT_Paras;         //parameter initial matrix
double *INIT_Metrics;        //initial results
double **Paras;              //intime parameters matrix
double *Metrics;             //intime metrics matrix
double *PARAS_Sum;           //sum of parameters of each simplex
double *HistMereics;         //metrics history

void downhill_simplex();
void read_config();
void read_init();
void update_status();
void get_paras_sum();
double update_simplex(int, double);
void check_paras_band(double *);
double get_metrics(double *);
void update_log(int, const char *);
int testfortermination();
void save_res();
void update_weight();


int
main(int argn, char **argv)
{
    int i;
    read_config();
	downhill_simplex();
	return 0;
}

void
read_config()
{
    int i, j;
    FILE *fp_init;
    time_t now;
    struct tm *timenow;
    char log_name[40];


    if(!(fp_init = fopen("downhill_config", "r"))){
        printf("downhill_config file does not exsit!\n");
        exit(1);
    }


    //read max iterator
    fscanf(fp_init, "%d", &NMAX);

    //read parameters dimension
    fscanf(fp_init, "%d", &NDIM);

    //read tolhist
    fscanf(fp_init, "%lf", &tolhist);

    //read benchmark configuration
    //fscanf(fp_init, "%s", &benchmark_path);
    fscanf(fp_init, "%s", &benchmark_case);
    fscanf(fp_init, "%s", &benchmark_res);
	fscanf(fp_init, "%s", &constraint_case);
	fscanf(fp_init, "%s", &constraint_res1);
	fscanf(fp_init, "%s", &constraint_res2);
    //printf("%s\n", benchmark_path);
    //printf("%s\n", benchmark_case);
    //printf("%s\n", benchmark_res);

    //allocate parameters range and intinial values
    RANGE_Paras  = (double **)malloc(NDIM * sizeof(double *));
    INIT_Paras   = (double **)malloc((NDIM + 1) * sizeof(double *));
    Paras        = (double **)malloc((NDIM + 1) * sizeof(double *));
    INIT_Metrics = (double *)malloc((NDIM + 1) * sizeof(double));
    Metrics      = (double *)malloc((NDIM + 1) * sizeof(double));
    PARAS_Sum    = (double *)malloc(NDIM * sizeof(double));
    HistMereics  = (double *)malloc(NMAX * sizeof(double));

    for(i = 0; i < NDIM; i++){
        RANGE_Paras[i] = (double *)malloc(2 * sizeof(double));
    }
    for(i = 0; i < NDIM + 1; i++){
        INIT_Paras[i] = (double *)malloc((NDIM) * sizeof(double));
        Paras[i]      = (double *)malloc((NDIM) * sizeof(double));
    } 

    fclose(fp_init);

	system("rm -rf simplex_log"); 
    //init log
#ifdef LOGFILE
    time(&now);
    timenow = localtime(&now);
    //sprintf(log_name, "simplex_log_%d.%d.%d-%d.%d.%d", timenow->tm_year+1900, \
            timenow->tm_mon, timenow->tm_mday, timenow->tm_hour, timenow->tm_min,\
            timenow->tm_sec);
	sprintf(log_name, "simplex_log");
    FP_Log  = fopen(log_name, "w");


    fprintf(FP_Log,"iterator+action+parameters+metrics\n");
    fflush(FP_Log);
    fclose(FP_Log);
#endif

    NUM_SIM = NDIM + 1; 
}


void 
downhill_simplex()
{	
    double rel_met;
    double ihi_save;
    int i, ihi, ilo, inhi, j, nfunk, min_all_shrunk;

    //init paras and metrics
    read_init();
    update_status();
	//begin iterate
	iterator = 0;
    while(testfortermination())
	{
        //read downhill configuration
        get_paras_sum();

		// find the highest, the next highest and the lowest points
        ilo  = 0;
        ihi  = 0;
        inhi = 0;
	    for(i = 0; i < NUM_SIM; i++){
            //printf("%e\n", Metrics[i]);
            if(Metrics[i] <= Metrics[ilo]){
                ilo = i;
            }
            if(Metrics[i] > Metrics[ihi]){
                ihi = i;
            }
        }
        for(i = 0; i < NUM_SIM; i++){
            if(Metrics[i] > Metrics[inhi] && i != ihi){
                inhi = i;
            }
        }

    

        //reflect the worst simplex, -1.0 = reflect
        rel_met = update_simplex(ihi, -1.0);
        //iterator++;
        update_log(-1, "reflect");

        //expand the reflect point
        if(rel_met <= Metrics[ilo]){
            rel_met = update_simplex(ihi, 2.0);
            HistMereics[iterator] = rel_met;
            iterator++;
            update_log(iterator, "expand");
        }
        //shrink the reflect point
        else if(rel_met >= Metrics[inhi]){
            ihi_save = Metrics[ihi];
            rel_met = update_simplex(ihi, 0.5);
            update_log(-1, "shrink");

            //Can’t seem to get rid of that high point. Better contract around the best point.    
            if(rel_met >= ihi_save){
                min_all_shrunk = (unsigned)(-1) >> 1; //INT MAX
                for(i = 0; i < NDIM + 1; i++){
                    if(i != ilo){
                        for(j = 0; j < NDIM; j++){
                            Paras[i][j] = 0.5 * (Paras[i][j] + Paras[ilo][j]);
                        }
                        check_paras_band(Paras[i]);
                        Metrics[i] = get_metrics(Paras[i]);
                        update_log(-1, "all_shrink");
                        if(Metrics[i] < min_all_shrunk){
                            min_all_shrunk = Metrics[i];
                        }
                   }
                }
                HistMereics[iterator] = min_all_shrunk;
                iterator++;
                update_log(iterator, "all_shrink");
                continue;
            }

            HistMereics[iterator] = rel_met;
            iterator++;
            update_log(iterator, "shrink");

        }
        //else keep the reflect point
        else{
            HistMereics[iterator] = rel_met;
            iterator++;
            update_log(iterator, "keep");
        }
	    update_weight();
    }
	
#ifdef LOGFILE
	//fclose(FP_Log);
#endif

    save_res();
}

void 
read_init()
{
    int i, j;
    FILE *fp_config;
    pw = 1.0;

    if(!(fp_config = fopen("subrange", "r"))){
        printf("%s: downhill_config file does not exsit!\n", __func__);
        exit(1);
    }

    //read parameters range and initial valuse
    for(i = 0; i < NDIM; i++){
        for(j = 0; j < 2; j++){
            fscanf(fp_config, "%lf", &RANGE_Paras[i][j]);
            //printf("%e ", RANGE_Paras[i][j]);
        }
        //printf("\n");
    }
    //printf("\n");
    for(i = 0; i < NDIM + 1; i++){
        for(j = 0; j < NDIM+1; j++){
			if(j != NDIM){
            	fscanf(fp_config, "%lf", &INIT_Paras[i][j]);
            	printf("%e ", INIT_Paras[i][j]);
           	}
			else{
				fscanf(fp_config, "%lf", &INIT_Metrics[i]);
			}
        }
        printf("\n");
    }

    for(i = 0; i < NDIM + 1; i++){
        //INIT_Metrics[i] = get_metrics(INIT_Paras[i]);
		printf("%e \n",INIT_Metrics[i]);
    }

    fclose(fp_config);
}

void
update_status(){
    //memcpy((void *)Paras, (void *)INIT_Paras, (NDIM + 1) * NDIM * sizeof(double));
    //memcpy((void *)Metrics, (void *)INIT_Metrics, (NDIM + 1) * sizeof(double));
    //printf("%e\n", Metrics[0]);
    int i, j;

    //reset parameters and metrics by configuration file
    for(i = 0; i < NDIM + 1; i++){
        for(j = 0; j < NDIM + 1; j++){
            if(j != NDIM){
                Paras[i][j] = INIT_Paras[i][j];
            }
            else{
                Metrics[i] = INIT_Metrics[i];
            }
        }
    }

    //reset update of configure file is 0
    //system("sed -i  '1s/1/0/g' subrange");
    //system("sed  '1s/1/0/g' downhill_config");
}

void
get_paras_sum()
{
    int i, j;

    memset(PARAS_Sum, 0, NDIM * sizeof(double));

    for(i = 0; i < NDIM; i++){
        for(j = 0; j < NDIM + 1; j++){
            PARAS_Sum[i] += Paras[j][i];
        }
    }
}

double 
update_simplex(int ihi, double fac)
{
    int i,j;
    double *paras_update;
    double fac1, fac2; 
    double metrics_update;

    paras_update = (double *)malloc(NDIM * sizeof(double));

    fac1 = (1.0 - fac)/NDIM;
    fac2 = fac1 - fac;

    for(i = 0; i < NDIM; i++){
        paras_update[i] = fabs(PARAS_Sum[i]*fac1 - Paras[ihi][i]*fac2);
    }

    check_paras_band(paras_update);
    metrics_update = get_metrics(paras_update);
   
    //If it’s better than the highest, then replace the highest.
    if(metrics_update < Metrics[ihi]){
        Metrics[ihi] = metrics_update;
        for(i = 0; i < NDIM; i++){
            PARAS_Sum[i] += paras_update[i] - Paras[ihi][i];
            Paras[ihi][i] = paras_update[i];
        }
    }

    free(paras_update);
    return metrics_update;
}

void 
check_paras_band(double* paras_array)
{
    int i;
    for(i = 0; i < NDIM; i++){
        if(paras_array[i] < RANGE_Paras[i][0]){
            paras_array[i] = RANGE_Paras[i][0];
        }
        else if(paras_array[i] > RANGE_Paras[i][1]){
            paras_array[i] = RANGE_Paras[i][1];
        }
    }
}

double 
get_metrics(double * paras_array)
{
    int i;
	double res;
    double bench_res;
    double cons_res1, cons_res2;
	FILE *fp_para_t;
    FILE *fp_res;
    FILE *fp_metrics;
    char bench_command[500];
    char bench_res_file[500];
    char cons_command[500];
    char cons_res_file[500];
	
	fp_para_t = fopen("sampling_data", "wr");
    fp_metrics = fopen("get_metrics", "aw");
	
    //sprintf(bench_command,  "%s", benchmark_case);
    //sprintf(bench_res_file, "%s", benchmark_res);
    //sprintf(cons_command,  "%s",  constraint_case);
    //sprintf(cons_res_file, "%s",  constraint_res1);

    // for bechmark
    //for(i = 0; i < NDIM; i++){
    //    sprintf(bench_command, "%s %e", bench_command, paras_array[i]);
    //    sprintf(cons_command, "%s %e", cons_command, paras_array[i]);
    //}

    //paras_array[0] = 0 - paras_array[0];
    for(i = 0; i < NDIM; i++){
        fprintf(fp_para_t, "%e ", paras_array[i]);
        fprintf(fp_metrics, "%e ", paras_array[i]);
    }
    fclose(fp_para_t); 

    printf("%s\n", benchmark_case);
    printf("%s\n", benchmark_res);
    printf("%s\n", constraint_case);
    printf("%s\n", constraint_res1);
    printf("%s\n", constraint_res2);

    system(benchmark_case);
    system(constraint_case);

    fp_res = fopen(benchmark_res, "r");
    fscanf(fp_res, "%lf", &bench_res);
    fclose(fp_res);

    fp_res = fopen(constraint_res1, "r");
    fscanf(fp_res, "%lf", &cons_res1);
    fclose(fp_res);
    
    fp_res = fopen(constraint_res2, "r");
    fscanf(fp_res, "%lf", &cons_res2);
    fclose(fp_res);

    printf("%lf\n", bench_res);
    printf("%lf\n", cons_res1);
    printf("%lf\n", cons_res2);

    for(i = 0; i < NDIM; i++){
        printf("%e ", paras_array[i]);
    }


    pw = 1000.0;
	if(cons_res1 > 1 &&  cons_res2 < 1){
		res = bench_res + fabs(cons_res1 - 1) * pw;
	}
    else if(cons_res2 > 1 &&  cons_res1 < 1){
		res = bench_res + fabs(cons_res2 - 1) * pw;
	}
    else if(cons_res1 > 1 &&  cons_res2 > 1){
		res = bench_res + fabs(cons_res1 - 1 + cons_res2 - 1) * pw;
	}
	else{
		res = bench_res;
	}
    printf(": %f : %f : %.4f \t### %f\n", bench_res, cons_res1, res, pw);
    fprintf(fp_metrics,"\t: %f : %f : %f : %f \t### %f\n", bench_res, cons_res1, cons_res2, res, pw);
    fclose(fp_metrics);
    //printf("%s: %s\n", command, __func__);
    //printf("%s: %s\n", res_file, __func__);
    return res;
}

void
update_log(int iterator, const char * action)
{
#ifdef LOGFILE
    int i, j;
    FP_Log  = fopen("simplex_log", "aw");
    fprintf(FP_Log, "%d \t %s\n", iterator, action);
    for(i = 0; i < NDIM + 1; i++){
        for(j = 0; j < NDIM + 1; j++){
            if(j != NDIM){
                fprintf(FP_Log, "%e ", Paras[i][j]);
            }
            else{
                fprintf(FP_Log, ": %e \n", Metrics[i]);
            }
        }
    }
    fflush(FP_Log); 
    fclose(FP_Log);
#endif
}

int
testfortermination()
{
    int i, n;
    int last_iterator;
    double last_max, last_min, tol_range;

    last_iterator = 10 + ceil(30. * NDIM / (NDIM + 1));

    //the last last_iterator is lower than tolhist
    if(iterator < NMAX && iterator > last_iterator){
        last_min = last_max = HistMereics[iterator - 1];
        for(i = iterator - 1; i >= (iterator - last_iterator - 1); i--){
            if(HistMereics[i] < last_min)
                last_min = HistMereics[i];
            if(HistMereics[i] > last_max)
                last_max = HistMereics[i];
        }
        tol_range = last_max - last_min;
        if(tol_range < tolhist)
            return 0; 
    }

    //iterator is greater than NMAX
    if(iterator > NMAX)
        return 0;

    //all of the simplex are equal
    n = 0;
    for(i = 0; i < NDIM; i++){
        if(Metrics[i] == Metrics[i+1]){
		   n++;
		}
    }

    if(n == NDIM)
        return 0;
    else
        return 1;
} 

void save_res()
{
    int i, j;
    FILE *fp_finial_res;
    
    fp_finial_res = fopen("final_res", "w");
    /*for(i = 0; i < NDIM; i++){
        fprintf(fp_finial_res, "%e ", Paras[0][i]);
    }*/
    fprintf(fp_finial_res, "%d\n", iterator);
    fprintf(fp_finial_res, "%e\n", Metrics[0]);
}

void update_weight()
{
	int i,j;
	double distance;
	double dnorm;
	double *paras_grav;
    
	paras_grav = (double *)malloc(NDIM * sizeof(double));

	//compute the center of gravity
	for(i=0; i < NDIM; i++){
		paras_grav[i] = PARAS_Sum[i] / (NDIM+1);
	}

    //compute the distance 
	distance = 0;
	for(i=0; i<NDIM+1; i++){
		for(j=0; j<NDIM; j++){
		   dnorm = sqrt((Paras[i][j] - paras_grav[j]) * (Paras[i][j] - paras_grav[j]) / (NDIM + 1));
		}
		distance += dnorm;
	}
	distance = distance / (NDIM + 1);

    //update the penalty
    if(pw < 1.0 / distance)
		pw = 1.0/distance; 	
    
}

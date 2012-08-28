#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include "R.h"

#define MAX_DIM 5000
#define MAX_DISTRIBUTION 100000


int main()
{
    return(0);
}


void DF(int dimy[1],float Xdata[MAX_DIM], float Ydata[MAX_DIM], int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1], int w_max[1], int w_max9[1], double df_sum_chi[1], double df9_sum_chi[1],
        double df_sum_like[1], double df9_sum_like[1], double df_max_chi[1], double df9_max_chi[1], double df_max_like[1], double df9_max_like[1]);
void Biometrika(int dim[1], double *X_distances, double *Y_distances,  int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1], int w_max[1], double sum_chi[1],
                 double sum_like[1], double max_chi[1], double max_like[1]);
void Create_Biometrika_Distribution(int dim[1], int montecarlo[1], double *X_distances, double *Y_distances,  int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1],
                                     int w_max[1], double *sum_chi, double *sum_like, double *max_chi, double *max_like, double *cumulative);
void Create_DF_Distribution(int dimy[1], int montecarloy[1], int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1], int w_max[1], int w_max9[1], double *v_df_sum_chi, double *v_df_sum_like,
                            double *v_df_max_chi, double *v_df_max_like, double *v_df9_sum_chi, double *v_df9_sum_like, double *v_df9_max_chi, double *v_df9_max_like, double *cumulative);
void Sequential(int dim[1], int montecarlo[1], double *X_distances, double *Y_distances,  int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1],
                                     int w_max[1], int score_id[1], double true_score[1], double A_threshold[1], double B_threshold[1], double p_val_threshold[1], double epsilon[1], double output_p_val[1], int output_monte[1]);


struct sample{
    int index;
    float data;
    };

int Inversions(int *permutation, int *source, int *inversion_count, int dim);
int Merge(int *permutation, int *source, int *inversion_count, int dim);

void Random_Permutation(int *permutation,int dim);
int compare (const void *a, const void *b);
int     compare_samples (const void *a, const void *b);
int  compare_samples (const void *a, const void *b)
     {
        struct sample a1=*(struct sample*)a;
        struct sample b1=*(struct sample*)b;
       return ((a1.data > b1.data) - (a1.data < b1.data));
     }

int compare (const void *a, const void *b)
{
    double a1 = *(const double*)a;
    double b1= *(const double*)b;
    if ( (a1 - b1) >0)
        return (1);
   else
        return(-1);
}
void Random_Permutation(int *permutation,int dim)
{
    int i,m, temp;
    for (i=0;i<dim;i++)
        permutation[i]=i;
    for (m=dim-1; m>0;m--){
        i=rand()%(m+1);
        temp=permutation[i];
        permutation[i]=permutation[m];
        permutation[m]=temp;
    }

}
void DF(int dimy[1],float Xdata[MAX_DIM], float Ydata[MAX_DIM], int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1], int w_max[1], int w_max9[1], double df_sum_chi[1], double df9_sum_chi[1],
        double df_sum_like[1], double df9_sum_like[1], double df_max_chi[1], double df9_max_chi[1], double df_max_like[1], double df9_max_like[1])
{
    int i,j,k,count[3][3],c1,c2,x[3],y[3], bad_sum, bad_max;
    double e, current_sum_chi, current_max_chi, current_like;
    double minx,miny,maxx,maxy;
    int dim=dimy[0];
    df_sum_chi[0]=0;
    df9_sum_chi[0]=0;
    df_sum_like[0]=0;
    df9_sum_like[0]=0;
    df_max_chi[0]=0;
    df9_max_chi[0]=0;
    df_max_like[0]=0;
    df9_max_like[0]=0;
    for (i=0;i<dim;i++)
        for (j=i+1;j<dim;j++){
            if (Xdata[i]<Xdata[j]){
                minx=Xdata[i];
                maxx=Xdata[j];
            }
            else {
                minx=Xdata[j];
                maxx=Xdata[i];
            }
            if (Ydata[i]<Ydata[j]){
                miny=Ydata[i];
                maxy=Ydata[j];
            }
            else {
                miny=Ydata[j];
                maxy=Ydata[i];
            }
            for(c1=0;c1<3;c1++)
                for (c2=0;c2<3;c2++)
                    count[c1][c2]=0;
            for (k=0;k<dim;k++)
                if ((k!=i) && (k!=j)){
                    if (Xdata[k]>maxx){
                            if (Ydata[k]>maxy)
                                count[2][2]++;
                            else{
                                if (Ydata[k]<miny)
                                    count[2][0]++;
                                else
                                    count[2][1]++;
                            }
                    }
                    else{
                        if (Xdata[k]<minx){
                            if (Ydata[k]>maxy)
                                count[0][2]++;
                            else{
                                if (Ydata[k]<miny)
                                    count[0][0]++;
                                else
                                    count[0][1]++;
                            }
                        }
                        else {//middle column
                        if (Ydata[k]>maxy)
                                count[1][2]++;
                            else{
                                if (Ydata[k]<miny)
                                    count[1][0]++;
                                else
                                    count[1][1]++;
                            }
                        }
                    }
                }
            for(c1=0;c1<3;c1++){
                x[c1]=0;
                y[c1]=0;
            }
            for (c1=0;c1<3;c1++){
                x[c1]=count[c1][0] + count[c1][1] + count[c1][2];
                y[c1]=count[0][c1] + count[1][c1] + count[2][c1];
            }
            //calculate df9
            bad_sum=0;
            bad_max=0;
            current_sum_chi=0;
            current_max_chi=0;
            current_like=0;
            for(c1=0;c1<3;c1++)
                for(c2=0;c2<3;c2++){
                    e=(double)(x[c1]*y[c2])/(double)(dim-2);
                    if (e>w_sum[0])
                        current_sum_chi+=(count[c1][c2]-e)*(count[c1][c2]-e)/e;
                    else
                        bad_sum=1;
                    if (e>w_max9[0])
                        current_max_chi+=(count[c1][c2]-e)*(count[c1][c2]-e)/e;
                    else
                        bad_max=1;

                    if (count[c1][c2]>0)
                        current_like+= count[c1][c2]*log((float)(count[c1][c2])/e);
                }
                if (bad_sum==0)
                    df9_sum_chi[0]+=current_sum_chi;
                if ( (current_max_chi > df9_max_chi[0]) && (bad_max==0) )
                    df9_max_chi[0]=current_max_chi;

                df9_sum_like[0]+=current_like;
                if (current_like > df9_max_like[0])
                    df9_max_like[0]=current_like;

                //calculate df
                e=(double)(x[1]*y[1])/(double)(dim-2);
                    if ( (e>w_sum[0]) && (dim-2-e>w_sum[0]))
                        df_sum_chi[0]+=(count[1][1]-e)*(count[1][1]-e)/e+(count[1][1]-e)*(count[1][1]-e)/(dim-2-e);
                    if ( (e>w_max[0]) && (dim-2-e>w_max[0]))
                        if (df_max_chi[0] < ((count[1][1]-e)*(count[1][1]-e)/e+(count[1][1]-e)*(count[1][1]-e)/(dim-2-e)))
                            df_max_chi[0]=(count[1][1]-e)*(count[1][1]-e)/e+(count[1][1]-e)*(count[1][1]-e)/(dim-2-e);

                    if ( (count[1][1] > 0) && (count[1][1]<dim-2)){
                        current_like= count[1][1]*log((float)(count[1][1])/e)+(dim-2-count[1][1])*log((float)(dim-2-count[1][1])/(dim-2-e));
                        df_sum_like[0]+= current_like;
                        if (current_like>df_max_like[0])
                            df_max_like[0]=current_like;
                    }

        }
}
void Biometrika(int dim[1], double *X_distances, double *Y_distances,  int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1], int w_max[1], double sum_chi[1], double sum_like[1], double max_chi[1],
                 double max_like[1])
{

    int dimy=dim[0], bad_sum, bad_max;
    double current_sum_chi, current_max_chi, current_like;
    double e;
    int i,j,k,c1,c2;
    int ***p_counts = malloc(dimy * sizeof (int **));
    if (p_counts == NULL){
        error("Memory TROUBLE29\n");
    }
    for (k=0;k<dimy;k++){
        p_counts[k]=malloc(2 * sizeof (int *));
        if (p_counts[k] == NULL){
            error("Memory TROUBLE30\n");
        }
        for(i=0;i<2;i++){
            p_counts[k][i]=malloc(2 * sizeof (int));
            if (p_counts[k][i]==NULL){
                error("Memory TROUBLE31\n");
            }
        }
    }

    struct sample *x_perm = (struct sample *) malloc(dimy * sizeof (struct sample));
    if (x_perm == NULL){
        error("TROUBLE33\n");
    }
    struct sample *y_perm = (struct sample *) malloc(dimy * sizeof (struct sample));
    if (y_perm == NULL){
        error("TROUBLE34\n");
    }
    int *y_rev = (int *) malloc((dimy-1) * sizeof (int));
    if (y_rev == NULL){
        error("TROUBLE35\n");
    }
    int *xy_perm = (int *) malloc((dimy-1) * sizeof (int));
    if (xy_perm == NULL){
        error("TROUBLE36\n");
    }
    int *inversion_count = (int *) malloc((dimy-1) * sizeof (int));
    if (inversion_count == NULL){
        error("TROUBLE37\n");
    }
    int *source = (int *) malloc((dimy-1) * sizeof (int));
    if (source == NULL){
        error("TROUBLE38\n");
    }
    int *xy_perm_temp = (int *) malloc((dimy-1) * sizeof (int));
    if (xy_perm_temp == NULL){
        error("TROUBLE39\n");
    }

    sum_chi[0]=0;
    max_chi[0]=0;
    sum_like[0]=0;
    max_like[0]=0;
    for (i=0;i<dimy;i++){
        for (j=0;j<dimy;j++){
            if (j<i){
                x_perm[j].index=j;
                x_perm[j].data=X_distances[i*dimy+j];
                y_perm[j].index=j;
                y_perm[j].data=Y_distances[i*dimy+j];
            }
            if (j>i){
                x_perm[j-1].index=j-1;
                x_perm[j-1].data=X_distances[i*dimy+j];
                y_perm[j-1].index=j-1;
                y_perm[j-1].data=Y_distances[i*dimy+j];
            }
        }
        qsort(x_perm,dimy-1, sizeof(x_perm[0]), compare_samples);
        qsort(y_perm,dimy-1, sizeof(x_perm[0]), compare_samples);
        for (j=0;j<dimy-1;j++){
            y_rev[y_perm[j].index]=j;
        }
        for (j=0;j<dimy-1;j++){
            xy_perm[j]=y_rev[x_perm[j].index];
            source[j]=j;
            inversion_count[j]=0;
            xy_perm_temp[j]=xy_perm[j];
        }
        Inversions(xy_perm_temp, source, inversion_count,dimy-1);
        for (j=0;j<dimy-1;j++){
            p_counts[j][1][0]=inversion_count[j];
            p_counts[j][1][1]=j-inversion_count[j];
            p_counts[j][0][1]=xy_perm[j]+inversion_count[j]-j;
            p_counts[j][0][0]=dimy-xy_perm[j]-inversion_count[j]-2;
        }
/*        for(j=0;j<dimy-1;j++){
            printf("%d ",p_counts[j][1][0]);
            printf("%d ",p_counts[j][1][1]);
            printf("%d ",p_counts[j][0][1]);
            printf("%d ",p_counts[j][0][0]);
            printf("\n");
        }
        printf("\n");
*/
        for (j=0;j<dimy-1;j++)
 /*           if (j!=i)*/{
/*                count[0][0]=0;
                count[0][1]=0;
                count[1][0]=0;
                count[1][1]=0;
                current_score=0;
                for(k=0;k<dimy;k++)
                    if ((j!=k) && (i!=k)){
                        count[(X_distances[i*dimy+j]>X_distances[i*dimy+k])][(Y_distances[i*dimy+j]>Y_distances[i*dimy+k])]++;
                }
                printf("%d %d %d %d\n",count[1][0], count[1][1], count[0][1], count[0][0]);
                */
                bad_sum=0;
                bad_max=0;
                current_sum_chi=0;
                current_max_chi=0;
                current_like=0;
                for(c1=0;c1<2;c1++)
                    for(c2=0;c2<2;c2++){
                        e=(double)(p_counts[j][c1][0]+p_counts[j][c1][1]) * (double)(p_counts[j][0][c2]+p_counts[j][1][c2])/(double)(dimy-2);
                    if (e>w_sum[0])
                        current_sum_chi+=(p_counts[j][c1][c2]-e)*(p_counts[j][c1][c2]-e)/e;
                    else
                        bad_sum=1;
                    if (e>w_max[0])
                        current_max_chi+=(p_counts[j][c1][c2]-e)*(p_counts[j][c1][c2]-e)/e;
                    else
                        bad_max=1;

                    if (p_counts[j][c1][c2]>0)
                        current_like+= p_counts[j][c1][c2]*log((float)(p_counts[j][c1][c2])/e);
                }
                if (bad_sum==0)
                    sum_chi[0]+=current_sum_chi;
                if ( (current_max_chi > max_chi[0]) && (bad_max==0) )
                    max_chi[0]=current_max_chi;

                sum_like[0]+=current_like;
                if (current_like > max_like[0])
                    max_like[0]=current_like;
            }
    }
 /* Freeing MEmory */
    for (k=0;k<dimy;k++){
        for(i=0;i<2;i++)
            free(p_counts[k][i]);
        free(p_counts[k]);
    }
    free(p_counts);
    free(x_perm);
    free(y_perm);
    free(y_rev);
    free(xy_perm);
    free(inversion_count);
    free(source);
    free(xy_perm_temp);


}
void Create_DF_Distribution(int dimy[1], int montecarloy[1], int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1], int w_max[1], int w_max9[1], double *v_df_sum_chi, double *v_df_sum_like,
                             double *v_df_max_chi, double *v_df_max_like, double *v_df9_sum_chi, double *v_df9_sum_like, double *v_df9_max_chi, double *v_df9_max_like, double *cumulative)
{
    int permutation[MAX_DIM],i,k;
    int dim=dimy[0];
    int monte_carlo_iterations=montecarloy[0];
    float Xdata[MAX_DIM], Ydata[MAX_DIM];
    double df_sum_chi[1],df9_sum_chi[1], df_sum_like[1], df9_sum_like[1], df_max_chi[1],df9_max_chi[1], df_max_like[1], df9_max_like[1];
    srand ( time(NULL) );

    for(i=0;i<monte_carlo_iterations;i++)
        cumulative[i]=(double)(i)/(double)(monte_carlo_iterations);
    for (i=0;i<dim;i++)
        Xdata[i]=i;
    for (k=0;k<monte_carlo_iterations;k++){
        Random_Permutation(permutation, dim);
        for(i=0;i<dim;i++)
            Ydata[i]=permutation[i];
        DF(dimy, Xdata, Ydata,sum_chi_flag, sum_like_flag, max_chi_flag, max_like_flag, w_sum, w_max, w_max9, df_sum_chi, df9_sum_chi, df_sum_like, df9_sum_like, df_max_chi, df9_max_chi, df_max_like, df9_max_like);
        if (sum_chi_flag[0]==1){
            v_df_sum_chi[k]=df_sum_chi[0];
            v_df9_sum_chi[k]=df9_sum_chi[0];
        }
        if (max_chi_flag[0]==1){
            v_df_max_chi[k]=df_max_chi[0];
            v_df9_max_chi[k]=df9_max_chi[0];
        }
        if (sum_like_flag[0]==1){
            v_df_sum_like[k]=df_sum_like[0];
            v_df9_sum_like[k]=df9_sum_like[0];
        }
        if (max_like_flag[0]==1){
            v_df_max_like[k]=df_max_like[0];
            v_df9_max_like[k]=df9_max_like[0];
        }
    }
    if (sum_chi_flag[0]==1){
        qsort(v_df_sum_chi, monte_carlo_iterations, sizeof (double), compare);
        qsort(v_df9_sum_chi, monte_carlo_iterations, sizeof (double), compare);
    }
    if (max_chi_flag[0]==1){
        qsort(v_df_max_chi, monte_carlo_iterations, sizeof (double), compare);
        qsort(v_df9_max_chi, monte_carlo_iterations, sizeof (double), compare);
    }
    if (sum_like_flag[0]==1){
        qsort(v_df_sum_like, monte_carlo_iterations, sizeof (double), compare);
        qsort(v_df9_sum_like, monte_carlo_iterations, sizeof (double), compare);
    }
    if (max_like_flag[0]==1){
        qsort(v_df_max_like, monte_carlo_iterations, sizeof (double), compare);
        qsort(v_df9_max_like, monte_carlo_iterations, sizeof (double), compare);
    }
}
void Sequential(int dim[1], int montecarlo[1], double *X_distances, double *Y_distances,  int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1],
                                     int w_max[1], int score_id[1], double true_score[1], double A_threshold[1], double B_threshold[1], double p_val_threshold[1], double epsilon[1], double output_p_val[1], int output_monte[1])
{
    int permutation[MAX_DIM],i,j,k,bad=1;
    int dimy=dim[0], good=1;
    int monte_carlo_iterationsy=montecarlo[0];
    double sum_chi[1], sum_like[1],max_chi[1],max_like[1], likelihood=0;
    double pa,p0, inc1, inc2, score=0;
    double *RY_distances = (double *) malloc(dimy*dimy * sizeof (double));
    if (RY_distances == NULL){
        error("Memory TROUBLE2\n");
    }
    pa=(1-epsilon[0])*p_val_threshold[0];
    p0=(1+epsilon[0])*p_val_threshold[0];
    inc1=log( (pa*(1-p0))/(p0*(1-pa)));
    inc2=log((1-pa)/(1-p0));
    srand ( time(NULL) );
    for (k=0; k < monte_carlo_iterationsy; k++ ){
        Random_Permutation(permutation, dimy);
        for (i=0; i<dimy; i++)
            for (j=0;j<dimy;j++)
                RY_distances[dimy*permutation[i]+permutation[j]]= Y_distances[dimy*i+j];
        if (score_id[0]<4)
            Biometrika(dim, X_distances, RY_distances,  sum_chi_flag, sum_like_flag, max_chi_flag, max_like_flag, w_sum, w_max, sum_chi, sum_like, max_chi, max_like);
        if (score_id[0]==0)
            score=sum_chi[0];
        if (score_id[0]==1)
            score=sum_like[0];
        if (score_id[0]==2)
            score=max_chi[0];
        if (score_id[0]==3)
            score=max_like[0];
        if (score_id[0]>3){
            error("Sequential suuports only Biometrika score at the moment\n");
        }
        if (score > true_score[0]){
            likelihood+=inc1;
            bad++;
        }
        likelihood+=inc2;
        if ( (likelihood<B_threshold[0]) && (good==1)){
            output_monte[0]=k+1;
            output_p_val[0]=1;
            break;
        }
        if (likelihood > A_threshold[0])
            good=0;
    }
    if (k==monte_carlo_iterationsy){
        output_monte[0]=k+1;
        output_p_val[0]=(double)(bad)/(double)(k+2);
    }
    free(RY_distances);

}
void Create_Biometrika_Distribution(int dim[1], int montecarlo[1], double *X_distances, double *Y_distances,  int sum_chi_flag[1], int sum_like_flag[1], int max_chi_flag[1], int max_like_flag[1], int w_sum[1],
                                     int w_max[1], double *v_sum_chi, double *v_sum_like, double *v_max_chi, double *v_max_like, double *cumulative)
{
    int permutation[MAX_DIM],i,j,k;
    int dimy=dim[0];
    int monte_carlo_iterationsy=montecarlo[0];
    double sum_chi[1], sum_like[1],max_chi[1],max_like[1];
    double *RY_distances = (double *) malloc(dimy*dimy * sizeof (double));
    if (RY_distances == NULL){
        error("Memory TROUBLE2\n");
    }
    srand ( time(NULL) );
    for (k=0; k < monte_carlo_iterationsy; k++ ){
        cumulative[k]=(double)(k)/(double)(monte_carlo_iterationsy);
        Random_Permutation(permutation, dimy);
        for (i=0; i<dimy; i++)
            for (j=0;j<dimy;j++)
                RY_distances[dimy*permutation[i]+permutation[j]]= Y_distances[dimy*i+j];
        Biometrika(dim, X_distances, RY_distances,  sum_chi_flag, sum_like_flag, max_chi_flag, max_like_flag, w_sum, w_max, sum_chi, sum_like, max_chi, max_like);
        if (sum_chi_flag[0]==1)
            v_sum_chi[k]=sum_chi[0];
        if (sum_like_flag[0]==1)
            v_sum_like[k]=sum_like[0];
        if (max_chi_flag[0]==1)
            v_max_chi[k]=max_chi[0];
        if (max_like_flag[0]==1)
            v_max_like[k]=max_like[0];
    }
        if (sum_chi_flag[0]==1)
            qsort(v_sum_chi, monte_carlo_iterationsy, sizeof (v_sum_chi[0]), compare);
        if (sum_like_flag[0]==1)
            qsort(v_sum_like, monte_carlo_iterationsy, sizeof (v_sum_like[0]), compare);
        if (max_chi_flag[0]==1)
            qsort(v_max_chi, monte_carlo_iterationsy, sizeof (v_max_chi[0]), compare);
        if (max_like_flag[0]==1)
            qsort(v_max_like, monte_carlo_iterationsy, sizeof (v_max_like[0]), compare);
    free(RY_distances);
}

int Inversions(int *permutation, int *source, int *inversion_count, int dim)
{
    if (dim==1)
        return 0;
    else{
        Inversions(permutation, source, inversion_count, dim/2);
        Inversions(&permutation[dim/2], &source[dim/2], inversion_count,dim-dim/2);
        Merge(permutation, source, inversion_count, dim);
    }
    return 0;
}

int Merge(int *permutation, int *source, int *inversion_count, int dim)
{
    int i;
    int left[MAX_DIM], right[MAX_DIM], left_source[MAX_DIM], right_source[MAX_DIM];
    int left_index=0, right_index=0;
    for (i=0;i<dim/2;i++){
        left[i]=permutation[i];
        left_source[i]=source[i];
    }
    for(i=0;i<dim-dim/2;i++){
        right[i]=permutation[i+dim/2];
        right_source[i]=source[i+dim/2];
    }
    for(i=0;i<dim;i++){
        if ( (left_index<dim/2) && (right_index<dim-dim/2)){
             if (left[left_index]<right[right_index]){
                permutation[i]=left[left_index];
                source[i]=left_source[left_index];
                left_index++;
            }
            else{
                permutation[i]=right[right_index];
                source[i]=right_source[right_index];
 //               printf("adding %d invs to %d\n", dim/2-left_index, source[i]);
                inversion_count[source[i]]+=(dim/2-left_index);
                right_index++;
            }
        }
        else{
            if (left_index<dim/2){
                permutation[i]=left[left_index];
                source[i]=left_source[left_index];
                left_index++;
            }
            if (right_index<dim-dim/2){
                permutation[i]=right[right_index];
                source[i]=right_source[right_index];
                right_index++;
            }

        }
    }
    return 0;
}

/*******************************************************************************************
 
 Copyright (c) 2020 Tyson LIttenberg & Neil Cornish
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
 **********************************************************************************************/

#include <hdf5.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_spline.h>

#include "Header.h"
#include "Constants.h"

double cosw(double t, double Tint);

#define DATASET "/sky/vgb/tdi"
#define FILENAME "LDC2_sangria_training_v2.h5"

static void forward_realft_wrapper(double *data, int n)
{
    /*
     Substitute forward GSL radix2 FFTs.
     
     E.g.,
    
     Replace: gsl_fft_real_radix2_transform(data, 1, n);
     With:    forward_realft_wrapper(data,n);
     
     Output arrays will be packaged as if
     the radix2 functions were used so that
     no other changes for a radix2-based
     code are needed.
     */
    
    gsl_fft_real_workspace * work = gsl_fft_real_workspace_alloc (n);
    gsl_fft_real_wavetable * real = gsl_fft_real_wavetable_alloc (n);
    double *scratch=malloc(n*sizeof(double));

    gsl_fft_real_transform(data, 1, n, real, work);

    // repack arrays as if we used radix2
    for(int i=0; i<n; i++)    scratch[i] = data[i];
    for(int i=1; i<n/2; i++)
    {
        data[i]   = scratch[2*i-1];
        data[n-i] = scratch[2*i];
    }

    free(scratch);
    gsl_fft_real_wavetable_free (real);
    gsl_fft_real_workspace_free (work);
}

struct TDI
{
    //Michelson
    double *X;
    double *Y;
    double *Z;
    
    //Noise-orthogonal
    double *A;
    double *E;
    double *T;
    
    //Number of data channels
    int Nchannel;
    
    //Number of frequency bins
    int N;
    
    //Data cadence
    double delta;
};


void free_tdi(struct TDI *tdi);
void alloc_tdi(struct TDI *tdi, int N, int Nchannel);

#define PI 3.141592653589793

// gcc -o segmentSangria segmentSangria.c -lhdf5 -lgsl

int main(int argc, char *argv[])
{
    int i, j, k, N, M, is, ii;
    int Ns, Ng;
    double fac;
    double t, f, x, y, Tend;
    char filename[1024];
    char command[1024];
    char strn[300];
    double alpha, s1, s2, Tcut;
    double *A, *E, *T, *Times;
    double *Arealft, *Erealft, *Trealft;
    double *AC, *EC;
    double *SA, *SE, *ST, *SN, *SNS;
    double SAE, SXYZ;
    
    FILE *in;
    FILE *out;
    
    struct TDI *tdi  = malloc(sizeof(struct TDI));
    
    /* LDC-formatted structure for compound HDF5 dataset */
    typedef struct tdi_dataset {
        double t;
        double    X;
        double    Y;
        double    Z;
    } tdi_dataset;
    static tdi_dataset *ldc_data;
    
    hid_t  file, dataset, dspace; /* identifiers */
    herr_t status; /* error handling */
    int ndims, Nrow, Ncol; /* dimension of dataset */
    double *data; /* array for dataset */
    
    /* Open an existing file. */
    file = H5Fopen(FILENAME, H5F_ACC_RDONLY, H5P_DEFAULT);
    
    /* Open an existing dataset. */
    dataset = H5Dopen(file, DATASET, H5P_DEFAULT);
    
    /* Get size of dataset */
    dspace = H5Dget_space(dataset);
    ndims = H5Sget_simple_extent_ndims(dspace);
    hsize_t dims[ndims];
    H5Sget_simple_extent_dims(dspace, dims, NULL);
    int Nsamples = dims[0];
    
    ldc_data = malloc(Nsamples*sizeof(struct tdi_dataset));
    alloc_tdi(tdi, Nsamples/2, 3);
    
    /* Setup memory for datatype handle */
    hid_t ldc_data_tid;
    
    ldc_data_tid = H5Tcreate(H5T_COMPOUND, sizeof(struct tdi_dataset));
    H5Tinsert(ldc_data_tid, "t", HOFFSET(struct tdi_dataset, t), H5T_IEEE_F64LE);
    H5Tinsert(ldc_data_tid, "X", HOFFSET(struct tdi_dataset, X), H5T_IEEE_F64LE);
    H5Tinsert(ldc_data_tid, "Y", HOFFSET(struct tdi_dataset, Y), H5T_IEEE_F64LE);
    H5Tinsert(ldc_data_tid, "Z", HOFFSET(struct tdi_dataset, Z), H5T_IEEE_F64LE);
    
    /* Read the dataset */
    status = H5Dread(dataset, ldc_data_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, ldc_data);
    
    /* Copy LDC-formatted structure into ldasoft TDI structure */
    for(i=0; i<Nsamples; i++)
    {
        double X = ldc_data[i].X;
        double Y = ldc_data[i].Y;
        double Z = ldc_data[i].Z;
        
        tdi->X[i] = X;
        tdi->Y[i] = Y;
        tdi->Z[i] = Z;
        
        //printf("%d %e\n", i, X);
        
        tdi->A[i] = (2.0*X-Y-Z)/3.0;
        tdi->E[i] = (Z-Y)/sqrt(3.0);
        tdi->T[i] = (X+Y+Z)/3.0;
        
    }
    
    
    tdi->delta = ldc_data[1].t - ldc_data[0].t;
    
    printf("Samples %d dt %f\n", Nsamples, tdi->delta);
    
    
    /* Close the dataset. */
    status = H5Dclose(dataset);
    
    /* Close the file. */
    status = H5Fclose(file);

    /* Free up memory */
    free(ldc_data);
    
    N = Nsamples;

    double Tseg, Tfull, Tdata;
    int Nseg, segs, seg;
    
    Tdata = (double)(Nsamples)*tdi->delta;
    Tfull = (double)(N)*tdi->delta;
    
    // default segment size is 30 days
    if(argc>1) Tseg = atof(argv[1])*day;
    else Tseg = 30.0*day;
    Nseg = (int)floor(Tseg/tdi->delta);
    segs = (int)floor(Tdata/Tseg);
    
    // allowing for 50% overlap
    segs = 2*segs-1;
    
    printf("Number of segments %d\n", segs);
    printf("Points per segement %d\n", Nseg);
    printf("Segment length (s) %.0f = %.0f (d)\n", Tseg, Tseg/day);
    
    Times = (double*)malloc(sizeof(double)*(N));
    
    /* FFT
    A = (double*)malloc(sizeof(double)*(N));
    E = (double*)malloc(sizeof(double)*(N));
    T = (double*)malloc(sizeof(double)*(N));
    */
    
    Arealft = (double*)malloc(sizeof(double)*(N));
    Erealft = (double*)malloc(sizeof(double)*(N));
    Trealft = (double*)malloc(sizeof(double)*(N));
    
    SN = (double*)malloc(sizeof(double)*(N));
    
    alpha = (2.0*t_tuke/Tseg);
    
    fac = sqrt(Tseg)/(double)(Nseg);   // Fourier scaling
    
    ii = 0;
    
    for (j = 0; j < 1; ++j)
    {
        
        for (i = 0; i < Nseg; ++i)
        {
            k = i+ii;
            
            /* FFT
            A[i] = tdi->A[k];
            E[i] = tdi->E[k];
            T[i] = tdi->T[k];
            */
            
            Arealft[i] = tdi->A[k];
            Erealft[i] = tdi->E[k];
            Trealft[i] = tdi->T[k];
        }
        

        sprintf(command, "AET_seg%d_t.dat", j);
        out = fopen(command,"w");
        for(i=0; i< Nseg; i++)
        {
            k = i+ii;
            t = (double)(k)*tdi->delta;
            //fprintf(out,"%.15e %.15e %.15e %.15e\n", t, A[i], E[i], T[i]);
            fprintf(out,"%.15e %.15e %.15e %.15e\n", t, Arealft[i], Erealft[i], Trealft[i]);
        }
        fclose(out);

        
        /* FFT
        tukey(A, alpha, Nseg);
        tukey(E, alpha, Nseg);
        tukey(T, alpha, Nseg);
         gsl_fft_real_radix2_transform(A, 1, Nseg);
         gsl_fft_real_radix2_transform(E, 1, Nseg);
         gsl_fft_real_radix2_transform(T, 1, Nseg);
        */
        
        tukey(Arealft, alpha, Nseg);
        tukey(Erealft, alpha, Nseg);
        tukey(Trealft, alpha, Nseg);

        forward_realft_wrapper(Arealft, Nseg);
        forward_realft_wrapper(Erealft, Nseg);
        forward_realft_wrapper(Trealft, Nseg);
        
        for(i=0; i< Nseg; i++)
        {
            /* FFT
            A[i] *= fac;
            E[i] *= fac;
            T[i] *= fac;
             */
            
            Arealft[i] *= fac;
            Erealft[i] *= fac;
            Trealft[i] *= fac;
        }
        
        /*
         prints Fourier domain data in a way that that can easily be plotted
         and compared
         */
       
        sprintf(command, "AET_seg%d_f_complex.dat", j);
        out = fopen(command,"w");
        for(i=1; i< Nseg/2; i++)
        {
            f = (double)(i)/Tseg;
            fprintf(out,"%e ", f);
            
            //radix2
            //fprintf(out,"%.16e %.16e %.16e %.16e %.16e %.16e\n", A[i], A[Nseg-i], E[i], E[Nseg-i], T[i], T[Nseg-i]);
            
            //realft
            fprintf(out,"%.16e %.16e %.16e %.16e %.16e %.16e\n", Arealft[i], Arealft[Nseg-i], Erealft[i], Erealft[Nseg-i], Trealft[i], Trealft[Nseg-i]);
        }
        fclose(out);
         /**/
        
        sprintf(command, "AET_seg%d_f.dat", j);
        out = fopen(command,"w");
        for(i=0; i< Nseg; i++)
        {
            f = (double)(i)/Tseg;
            
            //fprintf(out,"%e %.16e %.16e %.16e\n", f, A[i], E[i], T[i]);
            fprintf(out,"%e %.16e %.16e %.16e\n", f, Arealft[i], Erealft[i], Trealft[i]);

        }
        fclose(out);


        
        ii += Nseg/2;  //brick pattern
        
    }
    
    
    for (i = 0; i < N; ++i)
    {
        if(i < Nsamples)
        {
            /* FFT
            A[i] = tdi->A[i];
            E[i] = tdi->E[i];
            T[i] = tdi->T[i];
             */
            Arealft[i] = tdi->A[i];
            Erealft[i] = tdi->E[i];
            Trealft[i] = tdi->T[i];

        }
        else
        {
            /* FFT
            A[i] = 0.0;
            E[i] = 0.0;
            T[i] = 0.0;
             */
            Arealft[i] = 0.0;
            Erealft[i] = 0.0;
            Trealft[i] = 0.0;
        }
    }
    
    alpha = (2.0*t_tuke/Tdata);
    // we Tukey the original data length, not the padded
    
    /* FFT
    tukey(A, alpha, Nsamples);
    tukey(E, alpha, Nsamples);
    tukey(T, alpha, Nsamples);
     */
    
    tukey(Arealft, alpha, Nsamples);
    tukey(Erealft, alpha, Nsamples);
    tukey(Trealft, alpha, Nsamples);

    out = fopen("AET_t.dat","w");
    for (i = 0; i < N; ++i)
    {
        //fprintf(out,"%.15e %.15e %.15e %.15e\n", (double)(i)*tdi->delta, A[i], E[i], T[i]);
        fprintf(out,"%e %.16e %.16e %.16e\n", (double)(i)*tdi->delta, Arealft[i], Erealft[i], Trealft[i]);
        
    }
    fclose(out);
    
    /* FFT
    gsl_fft_real_radix2_transform(A, 1, N);
    gsl_fft_real_radix2_transform(E, 1, N);
    gsl_fft_real_radix2_transform(T, 1, N);
     */
    
    forward_realft_wrapper(Arealft, N);
    forward_realft_wrapper(Erealft, N);
    forward_realft_wrapper(Trealft, N);

    fac = sqrt(Tfull)/(double)(N);   // Fourier scaling
    
    for(i=0; i< N; i++)
    {
        /* FFT
        A[i] *= fac;
        E[i] *= fac;
        T[i] *= fac;
         */
        Arealft[i] *= fac;
        Erealft[i] *= fac;
        Trealft[i] *= fac;
    }
    
    out = fopen("AET_f.dat","w");
    
    for(i=0; i< N; i++)
    {
        f = (double)(i)/Tfull;
        
        //fprintf(out,"%e %.16e %.16e %.16e\n", f, A[i], E[i], T[i]);
        fprintf(out,"%e %.16e %.16e %.16e\n", f, Arealft[i], Erealft[i], Trealft[i]);
    }
    fclose(out);

    
    return 0;
    
}

void bwbpf(double *in, double *out, int fwrv, int M, int n, double s, double f1, double f2)
{
    /* Butterworth bandpass filter
     n = filter order 4,8,12,...
     s = sampling frequency
     f1 = upper half power frequency
     f2 = lower half power frequency  */
    
    if(n % 4){ printf("Order must be 4,8,12,16,...\n"); return;}
    
    int i, j;
    double a = cos(PI*(f1+f2)/s)/cos(PI*(f1-f2)/s);
    double a2 = a*a;
    double b = tan(PI*(f1-f2)/s);
    double b2 = b*b;
    double r;
    
    n = n/4;
    double *A = (double *)malloc(n*sizeof(double));
    double *d1 = (double *)malloc(n*sizeof(double));
    double *d2 = (double *)malloc(n*sizeof(double));
    double *d3 = (double *)malloc(n*sizeof(double));
    double *d4 = (double *)malloc(n*sizeof(double));
    double *w0 = (double *)malloc(n*sizeof(double));
    double *w1 = (double *)malloc(n*sizeof(double));
    double *w2 = (double *)malloc(n*sizeof(double));
    double *w3 = (double *)malloc(n*sizeof(double));
    double *w4 = (double *)malloc(n*sizeof(double));
    double x;
    
    for(i=0; i<n; ++i)
    {
        r = sin(PI*(2.0*(double)i+1.0)/(4.0*(double)n));
        s = b2 + 2.0*b*r + 1.0;
        A[i] = b2/s;
        d1[i] = 4.0*a*(1.0+b*r)/s;
        d2[i] = 2.0*(b2-2.0*a2-1.0)/s;
        d3[i] = 4.0*a*(1.0-b*r)/s;
        d4[i] = -(b2 - 2.0*b*r + 1.0)/s;
        w0[i] = 0.0;
        w1[i] = 0.0;
        w2[i] = 0.0;
        w3[i] = 0.0;
        w4[i] = 0.0;
    }
    
    for(j=0; j< M; ++j)
    {
        if(fwrv == 1) x = in[j];
        if(fwrv == -1) x = in[M-j-1];
        for(i=0; i<n; ++i)
        {
            w0[i] = d1[i]*w1[i] + d2[i]*w2[i]+ d3[i]*w3[i]+ d4[i]*w4[i] + x;
            x = A[i]*(w0[i] - 2.0*w2[i] + w4[i]);
            w4[i] = w3[i];
            w3[i] = w2[i];
            w2[i] = w1[i];
            w1[i] = w0[i];
        }
        if(fwrv == 1) out[j] = x;
        if(fwrv == -1) out[M-j-1] = x;
    }
    
    free(A);
    free(d1);
    free(d2);
    free(d3);
    free(d4);
    free(w0);
    free(w1);
    free(w2);
    free(w3);
    free(w4);
    
    return;
}

double tukeyt(double t, double tr, double Tint)
{
    double filter;
    
    filter = 1.0;
    if(t < tr) filter = 0.5*(1.0+cos(PI*(t/tr-1.0)));
    if(t > (Tint-tr)) filter = 0.5*(1.0+cos(PI*((Tint-t)/tr-1.0 )));
    if(t < 0.0) filter = 0.0;
    if(t > Tint) filter = 0.0;
    
    return filter;
    
}

double cosw(double t, double Tint)
{
    double filter;
    
    filter = sin(PI*(t/Tint));
    if(t < 0) filter = 0.0;
    if(t > Tint) filter = 0.0;
    
    return filter;
    
}


void tukey(double *data, double alpha, int N)
{
    int i, imin, imax;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    int Nwin = N-imax;
    
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i<imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(M_PI*( (double)(i-imax)/(double)(Nwin))));
        data[i] *= filter;
    }
    
}


void tukey_scale(double *s1, double *s2, double alpha, int N)
{
    int i, imin, imax;
    double x1, x2;
    double filter;
    
    imin = (int)(alpha*(double)(N-1)/2.0);
    imax = (int)((double)(N-1)*(1.0-alpha/2.0));
    
    int Nwin = N-imax;
    
    x1 = 0.0;
    x2 = 0.0;
    for(i=0; i< N; i++)
    {
        filter = 1.0;
        if(i<imin) filter = 0.5*(1.0+cos(M_PI*( (double)(i)/(double)(imin)-1.0 )));
        if(i>imax) filter = 0.5*(1.0+cos(M_PI*( (double)(i-imax)/(double)(Nwin))));
        x1 += filter;
        x2 += filter*filter;
    }
    x1 /= (double)(N);
    x2 /= (double)(N);
    
    *s1 = x1;
    *s2 = sqrt(x2);
    
}

void alloc_tdi(struct TDI *tdi, int NFT, int Nchannel)
{
    //Number of frequency bins (2*N samples)
    tdi->N = NFT;
    
    //Michelson
    tdi->X = calloc(2*tdi->N,sizeof(double));
    tdi->Y = calloc(2*tdi->N,sizeof(double));
    tdi->Z = calloc(2*tdi->N,sizeof(double));
    
    //Noise-orthogonal
    tdi->A = calloc(2*tdi->N,sizeof(double));
    tdi->E = calloc(2*tdi->N,sizeof(double));
    tdi->T = calloc(2*tdi->N,sizeof(double));
    
    int n;
    for(n=0; n<2*tdi->N; n++)
    {
        tdi->X[n] = 0.0;
        tdi->Y[n] = 0.0;
        tdi->Z[n] = 0.0;
        tdi->A[n] = 0.0;
        tdi->E[n] = 0.0;
        tdi->T[n] = 0.0;
    }
    
    //Number of TDI channels (X or A&E or maybe one day A,E,&T)
    tdi->Nchannel = Nchannel;
}

void free_tdi(struct TDI *tdi)
{
    free(tdi->X);
    free(tdi->Y);
    free(tdi->Z);
    free(tdi->A);
    free(tdi->E);
    free(tdi->T);
    
    free(tdi);
}


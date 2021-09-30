#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
//#include<athread.h>
#include"slave.h"
#include"dma.h"
#include"cpe_func.h"

void CopyEpiToFpWithCoordinateTransformCPE(void * args){
    int my_id = athread_get_id(-1);
    int n_tot_      = (int)(((unsigned long*)args)[0]);
    void * adr_epi_ = (void *)((unsigned long*)args)[1];
    void * adr_fp_  = (void *)((unsigned long*)args)[2];
    double PI = 3.14159265358979323846;
    //* Compute the task of each CPE
    int n_loc = n_tot_/NUMBER_OF_CPE + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot_/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot_ % NUMBER_OF_CPE) ? my_id : n_tot_ % NUMBER_OF_CPE );
    //* Memory allocation of local buffers
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_fp       = sizeof(fpLM);
    size_t bsize_epi       = sizeof(epiLM);
    size_t bsize_fp_array = bsize_fp * CHUNK_SIZE;
    size_t bsize_epi_array = bsize_epi * CHUNK_SIZE;
    epiLM * epi = (epiLM *) ldm_malloc( bsize_epi_array );
    fpLM *  fp =  (fpLM *)  ldm_malloc( bsize_fp_array );
    
    //* Compute moments on each CPE
    // (loop counters)
    int i,j,k; 
    // (local vars. for DMA comm.)
    dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    double one   = 1.0;
    double huge   = 1.0e300;

    /*
    double atanhi[] = {
        4.63647609000806093515e-01, // atan(0.5)hi 0x3FDDAC67, 0x0561BB4F
        7.85398163397448278999e-01, // atan(1.0)hi 0x3FE921FB, 0x54442D18
        9.82793723247329054082e-01, // atan(1.5)hi 0x3FEF730B, 0xD281F69B
        1.57079632679489655800e+00, // atan(inf)hi 0x3FF921FB, 0x54442D18
    };

    double atanlo[] = {
        2.26987774529616870924e-17, // atan(0.5)lo 0x3C7A2B7F, 0x222F65E2
        3.06161699786838301793e-17, // atan(1.0)lo 0x3C81A626, 0x33145C07
        1.39033110312309984516e-17, // atan(1.5)lo 0x3C700788, 0x7AF0CBBD
        6.12323399573676603587e-17, // atan(inf)lo 0x3C91A626, 0x33145C07
    };

    double aT[] = {
        3.33333333333329318027e-01, // 0x3FD55555, 0x5555550D
        -1.99999999998764832476e-01, // 0xBFC99999, 0x9998EBC4 
        1.42857142725034663711e-01, // 0x3FC24924, 0x920083FF 
        -1.11111104054623557880e-01, // 0xBFBC71C6, 0xFE231671
        9.09088713343650656196e-02, // 0x3FB745CD, 0xC54C206E 
        -7.69187620504482999495e-02, // 0xBFB3B0F2, 0xAF749A6D 
        6.66107313738753120669e-02, // 0x3FB10D66, 0xA0D03D51 
        -5.83357013379057348645e-02, // 0xBFADDE2D, 0x52DEFD9A 
        4.97687799461593236017e-02, // 0x3FA97B4B, 0x24760DEB 
        -3.65315727442169155270e-02, // 0xBFA2B444, 0x2C6A6C2F
        1.62858201153657823623e-02, // 0x3F90AD3A, 0xE322DA11 
    };
    */
    double * atanhi = (double *)ldm_malloc( sizeof(double)*4 );
    double * atanlo = (double *)ldm_malloc( sizeof(double)*4 );
    double * aT = (double *)ldm_malloc( sizeof(double)*11 );

    atanhi[0] = 4.63647609000806093515e-01; // atan(0.5)hi 0x3FDDAC67, 0x0561BB4F
    atanhi[1] = 7.85398163397448278999e-01; // atan(1.0)hi 0x3FE921FB, 0x54442D18
    atanhi[2] = 9.82793723247329054082e-01; // atan(1.5)hi 0x3FEF730B, 0xD281F69B
    atanhi[3] = 1.57079632679489655800e+00; // atan(inf)hi 0x3FF921FB, 0x54442D18

    atanlo[0] = 2.26987774529616870924e-17; // atan(0.5)lo 0x3C7A2B7F, 0x222F65E2
    atanlo[1] = 3.06161699786838301793e-17; // atan(1.0)lo 0x3C81A626, 0x33145C07
    atanlo[2] = 1.39033110312309984516e-17; // atan(1.5)lo 0x3C700788, 0x7AF0CBBD
    atanlo[3] = 6.12323399573676603587e-17; // atan(inf)lo 0x3C91A626, 0x33145C07

    aT[0] = 3.33333333333329318027e-01; // 0x3FD55555, 0x5555550D
    aT[1] = -1.99999999998764832476e-01; // 0xBFC99999, 0x9998EBC4 
    aT[2] = 1.42857142725034663711e-01; // 0x3FC24924, 0x920083FF 
    aT[3] = -1.11111104054623557880e-01; // 0xBFBC71C6, 0xFE231671
    aT[4] = 9.09088713343650656196e-02; // 0x3FB745CD, 0xC54C206E 
    aT[5] = -7.69187620504482999495e-02; // 0xBFB3B0F2, 0xAF749A6D 
    aT[6] = 6.66107313738753120669e-02; // 0x3FB10D66, 0xA0D03D51 
    aT[7] = -5.83357013379057348645e-02; // 0xBFADDE2D, 0x52DEFD9A 
    aT[8] = 4.97687799461593236017e-02; // 0x3FA97B4B, 0x24760DEB 
    aT[9] = -3.65315727442169155270e-02; // 0xBFA2B444, 0x2C6A6C2F
    aT[10] = 1.62858201153657823623e-02; // 0x3F90AD3A, 0xE322DA11 
    
    double y_x;
    double w,s1,s2,z_tmp;
    int ix,hx,id;
    
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_epi*nn);
        dma(dma_get, (long*)(((epiLM *)adr_epi_)+my_offset+i), (long*)(epi));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        for(j=0; j<nn; j++){
            const double x = epi[j].pos.x;
            const double y = epi[j].pos.y;
            const double z = epi[j].pos.z;

            double y_x = y / x;
            double phi = 0.0;

            //hx = __HI(x);
            hx = *(1+(int*)&y_x);
            ix = hx&0x7fffffff;
            if(ix>=0x44100000) {	// if |x| >= 2^66
                if(ix>0x7ff00000||
                   //(ix==0x7ff00000&&(__LO(x)!=0)))
                   (ix==0x7ff00000&&(*(int*)&y_x!=0)))
                    phi = y_x+y_x;		//
                if(hx>0) phi = atanhi[3]+atanlo[3];
                else     phi = -atanhi[3]-atanlo[3];
                id = -2;
            }
            else if (ix < 0x3fdc0000) {	// |x| < 0.4375
                if (ix < 0x3e200000) {	// |x| < 2^-29 
                    if(huge+x>one) phi = y_x;	// raise inexact
                }
                id = -1;
            }
            else {
                y_x = fabs(y_x);
                if (ix < 0x3ff30000) {		// |x| < 1.1875
                    //y_x = ((2.0*y_x-one)/(2.0+y_x))*(ix < 0x3fe60000) + ((y_x-one)/(y_x+one))*(ix >= 0x3fe60000);
                    //id = 0 + (ix >= 0x3fe60000);
                    if (ix < 0x3fe60000) {	// 7/16 <=|x|<11/16
                        id = 0;
                        y_x = (2.0*y_x-one)/(2.0+y_x); 
                    }
                    else {			// 11/16<=|x|< 19/16
                        id = 1;
                        y_x  = (y_x-one)/(y_x+one); 
                    }
                }
                else {
                    //y_x = ((y_x-1.5)/(one+1.5*y_x))*(ix < 0x40038000) + (-1.0/y_x)*(ix >= 0x40038000);
                    //id = 2*(ix < 0x40038000) + 3*(ix >= 0x40038000);
                    if (ix < 0x40038000) {	// |x| < 2.4375
                        id = 2;
                        y_x  = (y_x-1.5)/(one+1.5*y_x);
                    }
                    else {			// 2.4375 <= |x| < 2^66
                        id = 3;
                        y_x  = -1.0/y_x;
                    }
                }
            }
            /* end of argument reduction */
            z_tmp = y_x*y_x;
            w = z_tmp*z_tmp;
            /* break sum from i=0 to 10 aT[i]z**(i+1) into odd and even poly */
            s1 = z_tmp*(aT[0]+w*(aT[2]+w*(aT[4]+w*(aT[6]+w*(aT[8]+w*aT[10])))));
            s2 = w*(aT[1]+w*(aT[3]+w*(aT[5]+w*(aT[7]+w*aT[9]))));
            //s1 = z_tmp*(aT0+w*(aT2+w*(aT4+w*(aT6+w*(aT8+w*aT10)))));
            //s2 = w*(aT1+w*(aT3+w*(aT5+w*(aT7+w*aT9))));
            if (id == -1) phi = y_x - y_x*(s1+s2);
            else if(id >= 0){
                z_tmp = atanhi[id] - ((y_x*(s1+s2) - atanlo[id]) - y_x);
                phi = (hx<0)? -z_tmp:z_tmp;
            }
            if(x < 0.0){
                if(y >= 0.0){
                    phi = PI + phi;
                }
                else{
                    phi = phi - PI;
                }
            }
            
            //double phi = __ieee754_atan2(y, x);
            //double phi = atan(y/x);
            if(phi < 0.0) phi += 2.0*PI;
            else if(phi >= 2.0*PI) phi -= 2.0*PI;
            
            const double r = sqrt(x*x+y*y);
            fp[j].pos.x = phi;
            fp[j].pos.y = r;
            fp[j].pos.z = z;
            fp[j].vel.x = epi[j].vel.x;
            fp[j].vel.y = epi[j].vel.y;
            fp[j].vel.z = epi[j].vel.z;
            fp[j].id    = epi[j].id;
            fp[j].mass  = epi[j].mass;
        }

        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_fp*nn);
        dma(dma_put, (long*)(((fpLM *)adr_fp_)+my_offset+i), (long*)(fp));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    ldm_free(epi, bsize_epi_array);
    ldm_free(fp,  bsize_fp_array);    
}


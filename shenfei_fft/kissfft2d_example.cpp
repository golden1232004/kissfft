#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>
#include <unistd.h>
#include "kiss_fft.h"
#include "kiss_fftnd.h"

static int gaussian_shaped_labels(kiss_fft_cpx *pData,
                                  float *pTmp,
                                  float sigma, 
                                  int width, int height)
{
    if (0==pData)
        return -1;
    
    float *ptr=0;

    ptr = pTmp;
    float sigma2 = sigma*sigma*2;
    int cen_x = round(width/2.f)-1;
    int cen_y = round(height/2.f)-1;
    for (int h=0; h<height; h++)
    {
        int dy = h-cen_y;
        for (int w=0; w<width; w++)
        {
            int dx = w-cen_x;                
            *(ptr++) = (float)exp((-(dx*dx + dy*dy)/sigma2));
        }
    }
    
    // circshift the data    
    ptr = pTmp;
    for (int h=0; h<height; h++)
    {
        int idy = (h+cen_y)%height;
        ptr = pTmp + idy*width;
        for (int w=0; w<width; w++)  
        {    
            int idx = (w+cen_x)%(width);
            pData[0].r = ptr[idx];
            pData[0].i = 0;
            pData++;
        }
    }    
    return 0;
}

static
int getdims(int * dims, char * arg)
{
    char *s;
    int ndims=0;
    while ( (s=strtok( arg , ",") ) ) {
        dims[ndims++] = atoi(s);
        //printf("%s=%d\n",s,dims[ndims-1]);
        arg=NULL;
    }
    return ndims;
}

int main(int argc,char ** argv)
{
    if (3 != argc){
        printf("Input parameters:\n");
        printf(" eg-2d: -n 2,3\n");
        printf(" explaint: height=2, width=3\n");
        printf(" eg-3d: -n 2,3,4\n");
        printf(" explaint: channel=2, height=3, width=4");
        printf(" and data in memory should be order by channel.\n");
        return 0;
    }
    int k;
    int nfft[32];
    int ndims = 1;
    int isinverse=0;
    int numffts=1000,i;
    kiss_fft_cpx * buf;
    kiss_fft_cpx * bufout;
    kiss_fft_cpx * invbufout;
    int real = 0;

    nfft[0] = 1024;// default

    while (1) {
        int c = getopt (argc, argv, "n:ix:r");
        if (c == -1)
            break;
        switch (c) {
            case 'r':
                real = 1;
                break;
            case 'n':
                ndims = getdims(nfft, optarg );
                if (nfft[0] != kiss_fft_next_fast_size(nfft[0]) ) {
                    int ng = kiss_fft_next_fast_size(nfft[0]);
                    fprintf(stderr,"warning: %d might be a better choice for speed than %d\n",ng,nfft[0]);
                }
                break;
            case 'x':
                numffts = atoi (optarg);
                break;
            case 'i':
                isinverse = 1;
                break;
        }
    }
    int nbytes = sizeof(kiss_fft_cpx);
    for (k=0;k<ndims;++k)
        nbytes *= nfft[k];

    buf=(kiss_fft_cpx*)KISS_FFT_MALLOC(nbytes);
    bufout=(kiss_fft_cpx*)KISS_FFT_MALLOC(nbytes);
    invbufout = (kiss_fft_cpx*)KISS_FFT_MALLOC(nbytes);
    memset(buf,0,nbytes);
    float* pTmp = (float*)malloc(nfft[0]*nfft[1]*sizeof(float));
    float sigma = 1.5;
    int st = gaussian_shaped_labels(buf, pTmp, sigma, nfft[1], nfft[0]);
    if (3 == ndims){
        int i, j, k;
        for (i = 0; i < nfft[0]; i++){
            int ch = i * nfft[1]*nfft[2];
            for (j = 0; j < nfft[1]; j++){
                int row = ch + j*nfft[2];
                for (k = 0; k < nfft[2]; k++){
                    int idx = row + k;
                    buf[idx].r = idx;
                    buf[idx].i = 0;
                }
            }
        }
    }

    if (ndims==1) {
        if (real) {

        }else{
            kiss_fft_cfg st = kiss_fft_alloc( nfft[0] ,isinverse ,0,0);
            for (i=0;i<numffts;++i)
                kiss_fft( st ,buf,bufout );
            free(st);
        }
    }else{
        if (real) {
        }else{
            kiss_fftnd_cfg st= kiss_fftnd_alloc(nfft,ndims,isinverse ,0,0);
            kiss_fftnd( st ,buf,bufout );
            kiss_fftnd_cfg ist = kiss_fftnd_alloc(nfft, ndims, 1, 0, 0);
            kiss_fftnd(ist, bufout, invbufout);
            if (3 == ndims){
                int i, j;
                int val = nfft[0]*nfft[1]*nfft[2];
                for (i = 0; i < nfft[0]; i++){
                    int row = i * nfft[1];
                    for (j = 0; j < nfft[1]; j++){
                        int idx = row + j;
                        invbufout[idx].r /= val;
                        invbufout[idx].i /= val;
                    }
                }
            }
            free(st);
            free(ist);
        }
    }

    free(buf); free(bufout); free(invbufout);

    fprintf(stderr,"KISS\tnfft=");
    for (k=0;k<ndims;++k)
        fprintf(stderr, "%d,",nfft[k]);
    fprintf(stderr,"\tnumffts=%d\n" ,numffts);

    kiss_fft_cleanup();

    return 0;
}


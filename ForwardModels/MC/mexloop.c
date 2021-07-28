/*=================================================================
 * mexloop.c
 * 
 *  loop fuction in detPPL2hist.m 
 *
 * The calling syntax is:
 *
 *		[HIST] = mexloop(len, detID, bin, Track_atten)
 *      
 *      HIST: a vector [6000*1,
 *                      for detID = i,(for MATLAB index)
 *                      the bin: ((i-1)*60+1) : i*60  ]
 *
 *      len : length(Track_atten) [1*1 number of detected photons]
 *
 *      detID : a vector [len*1  detector ID integers frm 1 to 100 ]
 *       
 *      bin: a vector [len*1 time bin integers frm 1 to 60]
 *
 *      Track_atten: a vector [len*1 double]
 *      
 *      created by Jiang Jingjing 4/2/15.
 *
 *=================================================================*/

#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <mex.h>


/* Input Arguments */

#define	len_IN	prhs[0]
#define	detID_IN	prhs[1]
#define bin_IN  prhs[2]
#define Track_atten_IN  prhs[3]
/*#define num_det_IN  prhs[4] */
/* Output Arguments */

#define	hist_OUT	plhs[0]

#define N_DET 300
#define N_BIN 200

void mexhist(
             float	*hist,
             size_t len,
             float detID[],
             float	bin[],
             float Track_atten[]

             )
{
    size_t n_det, n_bin;
    int i;
    for(i=0; i<len; i++){
        n_det = (size_t)detID[i];
        n_bin = (size_t)bin[i];
        
        hist[(n_bin-1) * N_DET +  n_det -1] += Track_atten[i];
      
    }

    return;
}


void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    size_t len;

    float *HIST;
    float * detID;
    float * bin ;
    float *Track_atten;
 /*   size_t num_det;  */
    /* Check for proper number of arguments */
    
    if (nrhs != 4) {
	    mexErrMsgIdAndTxt( "MATLAB:mexloop:invalidNumInputs",
                "four input arguments required.");
    }
    /* input type control */
    if(mxGetClassID(detID_IN ) != mxSINGLE_CLASS)
        mexErrMsgTxt("input detID is not single, please use single(mat) as an input");
    if(mxGetClassID(bin_IN ) != mxSINGLE_CLASS)
        mexErrMsgTxt("input bin is not single, please use single(mat) as an input");
    if(mxGetClassID(Track_atten_IN ) != mxSINGLE_CLASS)
        mexErrMsgTxt("input Track_atten is not single, please use single(mat) as an input");
  
    /* Assign pointers to the various parameters */ 

    
    len = (int)(*mxGetPr(len_IN));
    detID = mxGetData(detID_IN);
    bin = mxGetData(bin_IN);
    Track_atten = mxGetData(Track_atten_IN);
  /*   num_det = (int)(*mxGetPr(num_det_IN));  */


  /* Create a matrix for the return argument */
     hist_OUT = mxCreateNumericMatrix((mwSize)N_DET ,(mwSize)N_BIN, mxSINGLE_CLASS, mxREAL); 
       HIST = (float*)mxGetPr(hist_OUT);

    /* Do the actual computations in a subroutine */
    mexhist(HIST,len ,detID,bin,Track_atten);
 
       return;
    
}

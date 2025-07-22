/*=================================================================
 * mexloop_uint32.c
 * 
 *  loop fuction create_hist.m
 *
 * The calling syntax is:
 *      
 *		HIST = mexloop_uint32(len, detID, bin, Track_atten, N_DET, N_BIN);   
 *   
 * Variables:
 *   len: scalar int (number of photons)
 *   detID: uint32 vector of detector IDs (1..N_DET)
 *   bin: uint32 vector of bin indices (1..N_BIN)
 *   Track_atten: single vector of photon intensities
 *   N_DET: scalar uint32 (total number of detectors)
 *   N_BIN: scalar uint32 (total number of bins)
 *   HIST: single matrix N_DET x N_BIN (detectors x bins)
 *
 * created by Jiang Jingjing 2025.07.22.
 *
 *=================================================================*/
#include "mex.h"
#include <string.h>  // for memset

/* Input Arguments */
#define LEN_IN          prhs[0]
#define DETID_IN        prhs[1]
#define BIN_IN          prhs[2]
#define TRACK_ATTEN_IN  prhs[3]
#define N_DET_IN        prhs[4]
#define N_BIN_IN        prhs[5]

/* Output Arguments */
#define HIST_OUT        plhs[0]

void mexhist(float *hist,
             size_t len,
             const uint32_t *detID,
             const uint32_t *bin,
             const float *Track_atten,
             size_t N_DET)
{
    size_t i;
    for (i = 0; i < len; i++) {
        size_t n_det = detID[i];
        size_t n_bin = bin[i];

        /* Validate indices inside loop can be skipped if MATLAB safeguards inputs,
           but adding here for safety */
        if (n_det == 0 || n_det > N_DET || n_bin == 0) {
            /* ignore out of range indices to avoid crash */
            continue;
        }

        /* indexing: column-major order.
         hist is N_DET x N_BIN matrix
         MATLAB element at (row = n_det, col = n_bin) is:
         hist[(n_det-1) + (n_bin-1)*N_DET]
        */
        hist[(n_det - 1) + (n_bin - 1) * N_DET] += Track_atten[i];
    }
}

void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Declare variables */
    size_t len;
    uint32_t *detID, *bin;
    float *Track_atten;
    size_t N_DET, N_BIN;
    float *hist;

    /* Check number of inputs */
    if (nrhs != 6) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:invalidNumInputs",
                          "Six input arguments required: len, detID, bin, Track_atten, N_DET, N_BIN.");
    }
    if (nlhs > 1) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:maxlhs",
                          "One output argument expected.");
    }

    /* Validate input types */
    if (!mxIsUint32(DETID_IN) || mxIsComplex(DETID_IN)) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:detIDType",
                          "detID must be a uint32 vector.");
    }
    if (!mxIsUint32(BIN_IN) || mxIsComplex(BIN_IN)) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:binType",
                          "bin must be a uint32 vector.");
    }
    if (!mxIsSingle(TRACK_ATTEN_IN) || mxIsComplex(TRACK_ATTEN_IN)) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:trackAttenType",
                          "Track_atten must be a single vector.");
    }
    if (!mxIsUint32(N_DET_IN) || mxIsComplex(N_DET_IN) || mxGetNumberOfElements(N_DET_IN) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:nDetType",
                          "N_DET must be a scalar uint32.");
    }
    if (!mxIsUint32(N_BIN_IN) || mxIsComplex(N_BIN_IN) || mxGetNumberOfElements(N_BIN_IN) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:nBinType",
                          "N_BIN must be a scalar uint32.");
    }
    if (!mxIsDouble(LEN_IN) || mxIsComplex(LEN_IN) || mxGetNumberOfElements(LEN_IN) != 1) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:lenType",
                          "len must be a scalar double.");
    }

    /* Get input data */
    len = (size_t)mxGetScalar(LEN_IN);

    /* Validate vector sizes */
    size_t detID_len = mxGetNumberOfElements(DETID_IN);
    size_t bin_len = mxGetNumberOfElements(BIN_IN);
    size_t track_len = mxGetNumberOfElements(TRACK_ATTEN_IN);

    if (len != detID_len || len != bin_len || len != track_len) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:lengthMismatch",
                          "len must match the length of detID, bin, and Track_atten vectors.");
    }

    detID = (uint32_t *)mxGetData(DETID_IN);
    bin = (uint32_t *)mxGetData(BIN_IN);
    Track_atten = (float *)mxGetData(TRACK_ATTEN_IN);

    N_DET = *((uint32_t *)mxGetData(N_DET_IN));
    N_BIN = *((uint32_t *)mxGetData(N_BIN_IN));

    if (N_DET == 0 || N_BIN == 0) {
        mexErrMsgIdAndTxt("MATLAB:mexloop_uint32:invalidDims",
                          "N_DET and N_BIN must be > 0.");
    }

    /* Create output matrix N_DET x N_BIN */
    HIST_OUT = mxCreateNumericMatrix((mwSize)N_DET, (mwSize)N_BIN, mxSINGLE_CLASS, mxREAL);
    hist = (float *)mxGetData(HIST_OUT);

    /* Initialize to zero */
    memset(hist, 0, sizeof(float) * N_DET * N_BIN);

    /* Compute histogram */
    mexhist(hist, len, detID, bin, Track_atten, N_DET);

    /* Return, hist_OUT assigned */

    return;
}

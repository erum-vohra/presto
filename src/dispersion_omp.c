#include "presto.h"

void dedisp_subbands(float *data, float *lastdata,
                     int numpts, int numchan,
                     int *delays, int numsubbands, float *result)
// De-disperse a stretch of data with numpts * numchan points into
// numsubbands subbands.  Each time point for each subband is a float
// in the result array.  The result array order is all the times for
// each subband, starting with lowest freq subband.  The delays (in
// bins) are in delays for each channel.  The input data and
// dispdelays are always in ascending frequency order.  Input data are
// contiguous channels with all of their points in time, starting with
// the lowest freq channel.
{
    const int chan_per_subband = numchan / numsubbands;
    long long ii, jj, kk, loffset;

    /* Initialize the result array */
    loffset = (long long)(numpts) * numsubbands;
    for (ii = 0; ii < loffset; ii++)
        result[ii] = 0.0;

    /* De-disperse into the subbands */
    for (ii = 0; ii < numchan; ii++) {
        const int subnum = ii / chan_per_subband;
        const int dind = delays[ii];
        float *sub = result + subnum * numpts;
        const long long loffset = ii * numpts;
        float *chan = lastdata + loffset + dind;

        for (jj = 0; jj < numpts - dind; jj++)
            sub[jj] += chan[jj];
        chan = data + ii * numpts;
        for (jj = numpts - dind, kk = 0; jj < numpts; jj++, kk++)
            sub[jj] += chan[kk];
    }
}


void float_dedisp(float *data, float *lastdata,
                  int numpts, int numchan,
                  int *delays, float approx_mean, float *result)
// De-disperse a stretch of data with numpts * numchan points. The
// delays (in bins) are in delays for each channel.  The result is
// returned in result.  The input data and delays are always in
// ascending frequency order.  Input data are ordered in time, with
// the channels stored together at each time point.
{
    long long ii, jj, kk;

    for (ii = 0; ii < numpts; ii++)
        result[ii] = -approx_mean;

    /* De-disperse */
    for (ii = 0; ii < numchan; ii++) {
        jj = ii + (long long)(delays[ii]) * numchan;
        for (kk = 0; kk < numpts - delays[ii]; kk++, jj += numchan)
            result[kk] += lastdata[jj];
        jj = ii;
        for (; kk < numpts; kk++, jj += numchan)
            result[kk] += data[jj];
    }
}


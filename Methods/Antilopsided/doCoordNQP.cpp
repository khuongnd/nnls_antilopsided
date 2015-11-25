#include <cmath>
#include <cstring>
#include <algorithm>
#include "pthread.h"
#include "math.h"
#include "mex.h"
#include "libs.h"
#include <time.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <ctime>

#define MAX_SUPPORTED_THREADS 100

using namespace std;
int FCDCycle(double **H, double *h, int k, double* df, double* x){
    For(iter, k){
        int pos = iter;
        /*int pos = -1, v;
        For(i, k)
            if ((df[i] < 0 || x[i] > 0) && H[i][i] > 0){
                dx = max(0.0, x[i] - df[i]/H[i][i]) - x[i];
                dfx = - 0.5 * dx *H[i][i] *  - df[i] * dx;
                if (maxValue < dfx){
                    pos = i;
                    maxValue = dfx;
                }
            }
        if (pos == -1) return iter; */
        if (H[pos][pos] > 0) {
            double dx = max(0.0, x[pos] - df[pos]/H[pos][pos]) - x[pos];
            x[pos] += dx;
            double *h = H[pos];
            For(i, k) df[i] += h[i] * dx;
        }
    }
    //mexPrintf("iter 2\n");
    return k;
}

void solveNQP(double *HH, double *h, double *x, int n, vector<double> &errors, 
        vector<double> &derivatives, vector<double>& times, int maxIter, double tolorance, double constant){
    //struct stableParams params;
    
    double **H = convert(HH, n, n);
    double *sqrtH = createV(n);
   
    double *df = createV(n);
    double *dbarf = createV(n);
    double *newX = createV(n);
    
    double den, a;
    double error, preDerivative, derivative;
    dotMV(H, x, n, n, df);
    addVV(df, h, n, df);
    removeKTTElements(df, x, n, dbarf);
    derivative = preDerivative = square(dbarf, n);
    
    error = (dotVV(df, x, n) + dotVV(h, x, n)) / 2 + constant;
    derivatives.push_back(derivative);
    errors.push_back(error);
    //times.push_back(getTime());
    
    For(iter, maxIter){
        
        //fast coordinate descent
        int nIter = FCDCycle(H, h, n, df, x);
        removeKTTElements(df, x, n, dbarf);
        derivative = square(dbarf, n);
        
        double error = (dotVV(df, x, n) + dotVV(h, x, n)) / 2 + constant;
        derivatives.push_back(derivative);
        errors.push_back(error);
        times.push_back(getTime());
        //mexPrintf("%d, %d, %.6e, %.6e, %.6e\n", iter, nIter, error, derivative, sum(x, n));
        
        //if (derivative > preDerivative)
        //    break;
        //preDerivative = derivative;
    }
    free(df);
    free(dbarf);
    free(newX);
    free(sqrtH);
}

///[xx, preDerivative, derivative, iter] = doIter (HH, hh, xx, maxIter, tolorance, verbose)
void mexFunction(int noOuts, mxArray *outs[], int noIns, const mxArray *inps[])
{
    vector<double> es, ds, ts;
    ts.push_back(getTime());
    
	double *H, *h;
	double tolorance, *iter, constant;
	int maxIter, n, verbose, a, b;

    //mexPrintf("%d, %d\n", noOuts, noIns);
	getArray(inps[0], &H, &n, &n);
    get1DArray(inps[1], &h, &a);
    if (a != n){
        //usage();
		mexPrintf("k=%d != k=%d, k=%d\n", a, n);
    }
    
	maxIter = (int)getDouble(inps[2]);
    tolorance = (double)getDouble(inps[3]);
    constant = (double)getDouble(inps[4]);
    
	/// Output arguments
	outs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	double* x =  mxGetPr(outs[0]);
    
    outs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    iter  = mxGetPr(outs[1]);
    
    solveNQP(H, h, x, n, es, ds, ts, maxIter, tolorance, constant);
    
    *iter = es.size();
    int it = (int)*iter;
    
    outs[2] = mxCreateDoubleMatrix(it, 1, mxREAL);
    double *errors  = mxGetPr(outs[2]);
    
    outs[3] = mxCreateDoubleMatrix(it, 1, mxREAL);
    double *ders  = mxGetPr(outs[3]);
    
    outs[4] = mxCreateDoubleMatrix(it, 1, mxREAL);
    double *times  = mxGetPr(outs[4]);
    
    For(i, it) errors[i] = es[i], ders[i] = ds[i], times[i] = ts[i];
}

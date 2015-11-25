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


void solveNQP(double *HH, double *h, double *x, int n, vector<double> &errors, 
        vector<double> &derivatives, vector<double>& times, int maxIter, double tolorance, double constant){
    //struct stableParams params;
    
    double **H = convert(HH, n, n);
    double *sqrtH = createV(n);
    double *saveDf = createV(n);
    double *saveX = createV(n);
    double *dx = createV(n);
    
    For(i, n) sqrtH[i] = (H[i][i] == H[i][i] && H[i][i] > 1e-9 ? 1.0/sqrt(H[i][i]) : 1.0);
    
    For(i, n) x[i] *= sqrtH[i], h[i] *= sqrtH[i];
    
    For(i, n) For(j, n) H[i][j] *= (sqrtH[i]*sqrtH[j]);
    
    double *df = createV(n);
    double *Qdx = createV(n), Qdf;
    double *dbarf = createV(n);
    double *newX = createV(n);
    
    double den, a;
    double error, preDerivative, derivative;
    dotMV(H, x, n, n, df);
    addVV(df, h, n, df);
    removeKTTElements(df, x, n, dbarf);
    derivative = preDerivative = square(dbarf, n);
    
    error = (dotVV(df, x, n) + dotVV(h, x, n)) / 2.0 + constant;
    derivatives.push_back(derivative);
    errors.push_back(error);
    
    For(iter, maxIter){
        assign(saveDf, df, n);
        assign(saveX, x, n);
        
        //antilopsided
        dotMV(H, dbarf, n, n, Qdx);
        den = dotVV(dbarf, Qdx, n);
        a = dotVV(dbarf, dbarf, n) / den;
        if (!(a != a || abs(a) < eps || abs(a) > 1e30)){
            For(i, n) {
                x[i] -= a * dbarf[i];
                df[i] -= a * Qdx[i];
                if (x[i] < 0) {
                    addkV(-x[i], H[i], n, df); 
                    x[i] = 0;
                }
            }
        }
        //fast coordinate descent
        FCDDf(H, h, n, df, x);
        
        //Accelerated Search
        For(i, n) dx[i] = saveX[i] - x[i];
        dotMV(H, dx, n, n, Qdx);
        den = dotVV(dx, Qdx, n);
        a = dotVV(df, dx, n) / den;
        if (!(a != a || abs(a) < eps || abs(a) > 1e30)){
            For(i, n) {
                x[i] -= a * dx[i];
                df[i] -= a * Qdx[i];
                if (x[i] < 0) {
                    addkV(-x[i], H[i], n, df); 
                    x[i] = 0;
                }
            }
        }
        
        //fast coordinate descent, choose max df[i]
        FCDDf(H, h, n, df, x);
        removeKTTElements(df, x, n, dbarf);
        derivative = square(dbarf, n);
        
        double error = (dotVV(df, x, n) + dotVV(h, x, n)) / 2.0 + constant;
        derivatives.push_back(derivative);
        errors.push_back(error);
        times.push_back(getTime(times[0]));
        
        preDerivative = derivative;
        printfFnc("iter = %d, time=%.4e, f(x)= %.5e,  df=%.5e\n", iter, getTime(times[0]), error, derivative);
    }
    For(i, n) x[i] *= sqrtH[i], h[i] *= sqrtH[i];
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
    
    mexPrintf("solve nqp\n");
    
    solveNQP(H, h, x, n, es, ds, ts, maxIter, tolorance, constant);
    
    *iter = es.size() - 1;
    int it = (int)*iter;
    
    outs[2] = mxCreateDoubleMatrix(it, 1, mxREAL);
    double *errors  = mxGetPr(outs[2]);
    
    outs[3] = mxCreateDoubleMatrix(it, 1, mxREAL);
    double *ders  = mxGetPr(outs[3]);
    
    outs[4] = mxCreateDoubleMatrix(it, 1, mxREAL);
    double *times  = mxGetPr(outs[4]);
    
    For(i, it) errors[i] = es[i+1], ders[i] = ds[i+1], times[i] = ts[i+1];
}

function res = FastNQP(XtX,Xty,opt, constant)
%FNNLS  Non-negative least-squares.
%
%       Adapted from NNLS of Mathworks, Inc.
%
%       x = fnnls(XtX,Xty) returns the vector X that solves x = pinv(XtX)*Xty
%       in a least squares sense, subject to x >= 0.
%       Differently stated it solves the problem min ||y - Xx|| if
%       XtX = X'*X and Xty = X'*y.
%
%       A default tolerance of TOL = MAX(SIZE(XtX)) * NORM(XtX,1) * EPS
%       is used for deciding when elements of x are less than zero.
%       This can be overridden with x = fnnls(XtX,Xty,TOL).
%
%       [x,w] = fnnls(XtX,Xty) also returns dual vector w where
%       w(i) < 0 where x(i) = 0 and w(i) = 0 where x(i) > 0.
%
%       See also NNLS and FNNLSb

%       L. Shure 5-8-87
%       Revised, 12-15-88,8-31-89 LS.
%       (Partly) Copyright (c) 1984-94 by The MathWorks, Inc.

%       Modified by R. Bro 5-7-96 according to
%       Bro R., de Jong S., Journal of Chemometrics, 1997, xx
%       Corresponds to the FNNLSa algorithm in the paper
%       http://newton.foodsci.kvl.dk/rasmus.html

%       Modified by S. Gunn 20-9-97
%
%  Reference:
%  Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.

% initialize variables
    tol = opt.eps
    if nargin < 3
        tol = 10*eps*norm(XtX,1)*length(XtX);
    end
    [m,n] = size(XtX);
    P = zeros(1,n);
    Z = 1:n;
    x = P';
    ZZ=Z;
    w = Xty-XtX*x;

    % set up iteration criterion
    iter = 0;
    itmax = 30*n;

    res.iter = 0;
    res.obj = [];
    res.d_barf = [];
    res.time = [];
    
    % outer loop to put variables into set to hold positive coefficients
    while any(Z) & any(w(ZZ) > tol) & toc(opt.startTime) < opt.maxTime
        res.iter = res.iter + 1;
        if (res.iter > opt.maxIter), break; end
        [wt,t] = max(w(ZZ));
        t = ZZ(t);
        P(1,t) = t;
        Z(t) = 0;
        PP = find(P);
        ZZ = find(Z);
        nzz = size(ZZ);
        z(PP')=(Xty(PP)'/XtX(PP,PP)');
        z(ZZ) = zeros(nzz(2),nzz(1))';
        z=z(:);
    % inner loop to remove elements from the positive set which no longer belong

        while any((z(PP) <= tol)) && iter < itmax

            iter = iter + 1;
            QQ = find((z <= tol) & P');
            alpha = min(x(QQ)./(x(QQ) - z(QQ)));
            x = x + alpha*(z - x);
            ij = find(abs(x) < tol & P' ~= 0);
            Z(ij)=ij';
            P(ij)=zeros(1,length(ij));
            PP = find(P);
            ZZ = find(Z);
            nzz = size(ZZ);
            z(PP)=(Xty(PP)'/XtX(PP,PP)');
            z(ZZ) = zeros(nzz(2),nzz(1));
            z=z(:);
        end
        x = z;
        w = Xty-XtX*x;
        df = XtX*x - Xty;
        BB = (x > 0 | df < 0);
        obj = x'*XtX*x/2.0 - Xty'*x + constant;
        d_barf = sum(df(BB).^2);
        
        res.obj = [res.obj obj];
        res.d_barf = [res.d_barf d_barf];
        res.time = [res.time toc(opt.startTime)];
        res.x = x;
        
        if (opt.verbose) 
            fprintf('%d,\t %E, \t%E\n',  res.iter, obj, d_barf);
        end
    end
end

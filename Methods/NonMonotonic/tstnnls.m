m=1e4;
n=1e3;
nargout = 2;
A = sprandn(m,n,.01);
xt = sprand(n,1,.2);
xt = full(xt);
b = A*xt;
opt=solopt;
opt.truex=1;

% for further speed uncomment the stuff below
opt.truex=0;
opt.verbose=0;
opt.xt=xt;
opt.x0 = 0.01*ones(n,1);
opt.maxTime = 16000;
size(A)
size(b)
size(x0)
out=bbnnls(A, b, opt); 
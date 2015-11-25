function opt = defaultOpt()
    
    opt.asgui = 0;
    opt.beta = 0.0498;
    opt.compute_obj = 1;
    % diminishing scalar; beta^0 =  opt.dimbeg
    % beta^k = opt.dimbeg / k^opt.dimexp
    opt.dimexp = .5;
    opt.dimbeg = 5;
    opt.maxit = 1000;
    opt.maxtime = 10;
    opt.maxnull = 10;
    opt.max_func_evals = 30;
    opt.pbb_gradient_norm = 1e-9;
    opt.sigma = 0.298;
    opt.step  = 1e-4;
    opt.tau = 1e-7;             
    opt.time_limit = 0;
    opt.tolg = 1e-3;
    opt.tolx = 1e-8;
    opt.tolo = 1e-5;
    opt.truex=0;
    opt.xt=[];
    opt.use_kkt = 0;
    opt.use_tolg = 1;
    opt.use_tolo = 0;
    opt.use_tolx = 0;
    opt.useTwo = 0;
    opt.verbose = 1;                    % initially
    if nargin == 1
      opt.variant = varargin{1};
    else   % Default
      opt.variant = 'SBB';
    end

    opt.maxIter = 1000;
    opt.verbose = 1;
    opt.eps = 1e-10;
    opt.tolorance = 1e-5;
end
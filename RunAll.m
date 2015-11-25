addpath (genpath('Methods'));

files = GetFiles();

funcs = { @AloCoordNNLS, @CoordNNLS, @AcceleratedNNLS, @FastNNLS, @bbnnls, @FCD};

maxNumCompThreads(1)

for j=1:length(funcs),
    errors = [];
    for i=5:6,
        file = sprintf('datasets/data_%d.mat', i);
        data = load(file);
        
        A = data.A;
        b = data.b;
        [d, n] = size(A);
        x0 = zeros([n, 1]);
        
        opt = defaultOpt();
        opt.maxIter = 3000;
        opt.eps = 1e-10;
        opt.verbose = 0;
        
        opt.xt = x0;
        opt.maxit = opt.maxIter;
        opt.tolg = opt.eps;
        opt.x0 = x0;
        opt.maxTime = 1600;
        opt.verbose = 0;
        opt.accuracy = 1;

        from = strfind(files{i}, '/');
        to = strfind(files{i}, '.');
        
        res = funcs{j}(A, b, opt);
        errors = [errors res.finalObj];
        sum(res.x > 0)
        
        fileOut = sprintf('Results/data_%d_%s.log', i, func2str(funcs{j}));
        
        fprintf('%s\n', fileOut);
        
        writeMatrix(fileOut, full([res.time' res.obj' res.d_barf']), ',');
    end
    fileOut = sprintf('Results/errors_%s.log', func2str(funcs{j}));
    writeMatrix(fileOut, errors', ',');
end

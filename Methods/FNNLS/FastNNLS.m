function res = FastNNLS(A,b,opt)
    startTime = tic;
    opt.startTime = startTime;
    
    res = FastNQP(A'*A, A'*b, opt, b'*b/2.0);
    optimal = sum((A*res.x - b).^2)/2.0;
    
    res.obj = [res.obj optimal];
    res.d_barf = [res.d_barf toc(startTime)];
    res.time = [res.time toc(startTime)];
    
    res.finalObj = norm(A*res.x - b, 'fro')^2/2;
    
    if (opt.verbose) 
        fprintf('Optimal: %.20E, \t, %.20E\n', optimal, optimal - res.obj(length(res.obj)-1));
    end
end
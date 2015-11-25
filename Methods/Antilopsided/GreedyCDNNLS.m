function res = GreedyCDNNLS(A,b,opt)
    startTime = tic;
    Q = full(A'*A);
    q = full(-b'*A)';
    
    constant = b'*b / 2;
    initTime = toc(startTime);
    [x, iter, errors, ders, times] = doGreedyCDNQP(Q, q, opt.maxIter, opt.tolorance, constant);
    
    res.iter = iter;
    res.obj = errors';
    res.d_barf = ders';
    res.time = (times - min(times) + initTime)';
    res.finalObj = norm(A*x - b, 'fro')^2/2.0;
    res.x = x;
    
    for i=1:iter,
         fprintf('%d,\t %E\t, %E\t, %E\n', i, res.time(i), res.obj(i), res.d_barf(i));
    end
    fprintf('Optimal: %.20E\n', norm(A*x - b, 'fro')^2/2.0);
end
function res = AcceleratedNNLS(A, b, opt)
    %    min    0.5*||Ax-b||^2, s.t. x >= 0
    startTime = tic;

    Q = A'*A;
    q = (-b'*A)';
    constant = b'*b / 2;
    %clear H;
    % the NNLS is equivalent to min 1/2*x'*Q*x+q'*x

    res.iter = 0;
    res.d_barf = [];
    res.obj = [];
    if isfield(opt, 'x0')
        x = opt.x0;
    else
        x = zeros(size(Q, 1), 1);
    end
    df = Q*x + q;
    res.obj = [];
    res.d_barf = [];
    res.time = [];

    L = sqrt(sum(sum(Q.^2)));
    Z=x;
    H=Z;    % Initialization
    Grad=Q*Z-q;     % Gradient
    alpha1=1;
    iter = 0;
    for iter = 1:opt.maxIter,
        if toc(startTime) > opt.maxTime,
            break;
        end;
        
        H0=H;
        H=max(Z-Grad/L,0);    % Calculate sequence 'Y'
        alpha2=0.5*(1+sqrt(1+4*alpha1^2));
        Z= max(0, H + ((alpha1-1)/alpha2) * (H-H0));
        alpha1=alpha2;

        obj = H'*Q*H/2.0 + q'*H + constant;
        Grad = Q*H+q;

        res.x = H;
        BB = (H > 0 | Grad < 0);

        d_barf = Grad(BB);
        
        if opt.accuracy,
            obj = sum((A*res.x - b).^2/2.0);
        end

        res.obj = [res.obj obj];
        res.d_barf = [res.d_barf d_barf'*d_barf];
        res.time = [res.time toc(startTime)];

        if (opt.verbose)
            fprintf('%d,\t, %E\t, %E\n', iter, obj, d_barf'*d_barf);
        end
        if (d_barf'*d_barf < opt.eps)
            break;
        end

        Grad=Q*Z+q;
    end

    
    res.d_barf = [res.d_barf toc(startTime)];
    res.obj = [res.obj sum((A*res.x - b).^2)/2.0];
    res.finalObj = norm(A*res.x - b, 'fro')^2/2.0;
    res.time = [res.time toc(startTime)];
    res.finalObj = norm(A*x - b, 'fro')^2/2;
    return;
end
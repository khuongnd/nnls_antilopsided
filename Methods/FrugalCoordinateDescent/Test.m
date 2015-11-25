addpath (genpath('..'));
rng(3579);
for i=1:1,
    switch i
        case 1,
            d = 80;
            n = 100;
            A = rand(d, n) * 10.0;
            for k=1:n/4,
                i1 = min(int32(rand()*n)+1, n);
                i2 = min(int32(rand()*n)+1, n);
                v1 = A(:, i1);
                v2 = rand(d, 1) * 2.0;
                v2 = v1 + v2;
                A(:, i2) = v2;
            end
        case 2,
            d = 15000;
            n = 10000;
            A = rand(d, n) * 10.0;
            for k=1:n/4,
                i1 = min(int32(rand()*n)+1, n);
                i2 = min(int32(rand()*n)+1, n);
                v1 = A(:, i1);
                v2 = rand(d, 1) * 2.0;
                v2 = v1 + v2;
                A(:, i2) = v2;
            end
        case 3,
            data = load('../../datasets/cifar.mat');
            A = double(data.data);
        case 4,
            data = load('../../datasets/20News.mat');
            A = data.fea(1:10000,:)'*10;
        case 5,
            d = 20000;
            n = 10000;
            A = sprand(d, n, 0.1) * 10.0;
            for k=1:n/4,
                i1 = min(int32(rand()*n)+1, n);
                i2 = min(int32(rand()*n)+1, n);
                v1 = A(:, i1);
                v2 = rand(d, 1) * 10;
                v2(v1 == 0) = 0;
                v2 = sparse(v2);
                v2 = v1 + sparse(v2);
                A(:, i2) = sparse(v2);
            end
            A = sparse(A);
        case 6,
            data = load('../../datasets/cifar.mat')';
            data = data.data(1:10000,:)';
            A = double(data);
    end;
    [d, n] = size(A);
    x = rand(n, 1);
    x (x < 0.2) = 0;
    b = A*x*10*max(max(A)) + rand(d, 1) * 100;
    size(b)
    save(sprintf('data_%d.mat', i), 'A', 'b');
    size(A)
    
    clear opt;
    opt = defaultOpt();
    opt.maxIter = 200;
    opt.eps = 1e-10;
    opt.tolorance = 1e-4;
    opt.verbose = 0;

    x0 = zeros(size(A,2), 1);
    opt.xt=x0;
    opt.tolg = opt.eps;
    opt.x0 = x0;
    opt.maxTime = 1600;
    opt.verbose = 0;
    opt.accuracy = 1;

    %maxNumCompThreads(4)

    his = FCD(A, b, opt);
    %res = CoordNNLS(A, b, opt);
    
    fprintf('%f, %f, %f\n', his.time, his.d_barf,  his.obj);
    size(his.time)
    size(his.d_barf)
    size(his.obj)
    %res.time'
end
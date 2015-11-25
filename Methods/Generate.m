function [A, b] = Generate(n, m, sparse, normalized, l)
    A = double(zeros(n, m));
    for i=1:m,
        A(:, i) = l*normalize(generateSparse(n, sparse), normalized);
    end
    x = l*generateSparse(m, sparse);
    b = A*x + rand(n, 1);
end

function x = generateSparse(n, sparse)
    x = rand(n, 1);
    for i=1:1:(n*sparse),
        x(i, 1) = 0.0;
    end
end

function res = normalize(x, normalized)
    if normalized == 1,
        x = x / norm(x, 'fro');
    end
    res = x;
end
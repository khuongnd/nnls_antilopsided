function [A b] = ReadInput(file)
A = ReadMatrix.readMatrix(file);
[n dimension] = size(A);
b = A(n,:)';
A = A(1:n-1,:)';

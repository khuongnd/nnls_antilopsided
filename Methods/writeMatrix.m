function result = writeMatrix(filename, A, seperator)
    [M, N] = size(A);
    fileID = fopen(filename,'W');
    for i=1:M
        for j=1:N
            if (j < N)
                fprintf(fileID,'%.20e%s', A(i,j), seperator);
            else
                fprintf(fileID,'%.20e\n', A(i,j));
            end
        end
    end
    result = true;
    fclose(fileID);
end
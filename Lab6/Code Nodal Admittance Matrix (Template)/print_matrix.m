function print_matrix(M,name,precision)


[n_rows,n_columns] = size(M);

% diagonal elements

for i=1:min(n_rows,n_columns)
    disp([name '(' int2str(i) ',' int2str(i) '): ' num2str(M(i,i),precision)]);
end

% off-diagonal elements

for i=1:n_rows
    for j=(i+1):n_columns
        if(abs(M(i,j))>1e-6)
            disp([name '(' int2str(i) ',' int2str(j) '): ' num2str(M(i,j),precision)]);
        end
    end
end

end
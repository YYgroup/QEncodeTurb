function mat_out = rowExpandMat(mat_in,K)
%Inputs an m x n matrix and an integer K, where mod(n,K) must equal 0, and
%expands each of the m rows K times. K is assumed here to represent the
%dimensionality of the leading index of each of the individual m rows.

[m,n] = size(mat_in);
mat_out = zeros(m*K,n/K);

for i = 1:m
    j = 1 + (i-1)*K;
    mat_out(j:(j-1+K),:) = reshape(mat_in(i,:),[K,n/K]);
end

end

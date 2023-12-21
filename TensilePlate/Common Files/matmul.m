function C = matmul(A,B,m,n,p)
% multiply of two tensor A_{m,n} * B_{n,p}
% A: m*n
% B: n*p
% C: m*p
% It is the same as A_{m,n} * B_{n,p} A_,B_ is the matrix form of A,B
np = size(A,1);
C = zeros(np,m*p);
for i = 1:m
    for j = 1:n
        for k = 1:p
            i1 = (i-1)*n + j;
            i2 = (j-1)*p + k;
            i3 = (i-1)*p + k;
            C(:,i3) = C(:,i3)+A(:,i1).*B(:,i2);
        end
    end
end
end
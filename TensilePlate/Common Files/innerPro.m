function K = innerPro(wf, X, Y, vj, ih,totnode,dof)
% K - dimension: (totnode, dof^2)
% obtain the inner product of two states X and Y
% k = (11 12)
%     (21 22)
% K = (11,12,21,22)
% Tensor product,it is the same as Y(i,:)'*X(i,:)
m = size(X,2);
n = size(Y,2);
K = zeros(totnode,dof^2);
for i = 1:m
    for j = 1:n
        k = (i-1)*n + j;
        K(:,k) = accumarray(ih, wf.*X(:,i).*Y(:,j).*vj);
    end
end
end
function [A,B] = spectral_graph_partitioning(L,balanced)
%L symetric sparse GGL matrix
d = diag(L); %degree
n = length(d);
P = spdiags(sqrt(1./d),0,n,n);

W2 = -L; %copy L into W2
W2 = spdiags(zeros(n,1),0,W2); %make diagonal of W2 equal to zero

W2 = (P*W2*P); %normalize
W2 = (W2 + W2')/2; %make sure symmetric
D2 = spdiags(sum(W2)',0,n,n);

[evect, ~] = eigs(D2-W2,1);
if (~balanced)
    B = (evect<0);
else
    sorted_evect = sort(evect);
    th = sorted_evect(round(n/2));
    B = (evect<=th);
end
%now make sure that set A has more nodes than B
if(nnz(B)>n/2)
    A=B;
    B = ~A;
else
    A = ~B;
end
end
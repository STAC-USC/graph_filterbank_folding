% Authors - Eduardo Pavez <eduardo.pavez.carvelli@gmail.com, pavezcar@usc.edu>
% Copyright Eduardo Pavez,  University of Southern California
%Los Angeles, USA, 02/11/2023
% E. Pavez, B. Girault, A. Ortega, and P. A. Chou. 
%"Two channel filter banks on arbitrary graphs with positive semi-definite variation operators". 
%IEEE Transactions on Signal Processing 2023
%https://arxiv.org/abs/2203.02858


%line/path graph example
%this code displays the (M,Q)-GFT for a line graph (bipartite) and a
%non-bipartite line graph
clear;
n = 11;%number of nodes of line graph

%first construct line/path graph (bipartite)
W = zeros(n,n);
for i=1:n-1
    W(i,i+1)=1;
end
W = spones(W + W');
D = diag(sum(W));
L = D-W;
%Q = D;%because graph is bipartite Q=D
Z = D\L;
%B = speye(n) - Z;

%now add non bipartite edges (see figure in paper)

W2 = W;

W2(1,3)=1;
W2(2,4)=1;
W2(5,7)=1;
W2(6,8)=1;
W2(9,11)=1;

W2 = spones(W2+W2');%symmetrize and make sure all edges have unit weight

D2 = diag(sum(W2));
L2 = D2-W2;
Q2 = L2;
A = [1 3 5 7 9 11];%odd nodes
B = [2 4 6 8 10];%even nodes
Q2(A,B)=0;
Q2(B,A)=0;
Z2 = Q2\L2;

%B2 = speye(n) - Z2;

% now compute and plot the GFT


[Ud,Vd] = eigs(L,D,n,'smallestabs'); %bipartite line graph (L,D)-GFT (with spectral folding)
[U2d,V2d] = eigs(L2,D2,n,'smallestabs'); %non bipartite line graph (L2,D2)-GFT (without spectral folding)
[U2q,V2q] = eigs(L2,Q2,n,'smallestabs'); %non bipartite line graph (L2,Q2)-GFT (with spectral folding)

lambdas = lambda_title(n);

plot_MQ_GFT(Ud,Vd,n,lambdas);
saveas(gcf,'line_bipartite_LD.png','png')
figure;
plot_MQ_GFT(U2d,V2d,n,lambdas);
saveas(gcf,'line_nonbipartite_LD.png','png')
figure;
plot_MQ_GFT(U2q,V2q,n,lambdas);
saveas(gcf,'line_nonbipartite_LQ.png','png')
function [] = plot_MQ_GFT(U,V,N,lambdas)
lambda=diag(V);
lambda(1)=0;
m = floor(N/2);%N is odd

for p=1:m
    subplot(2,m+1,p);
    u = U(:,p);
    if(u(1)^2>0)
        u = u.*sign(u(1));
    end
    stem(u,'b','linewidth',1);
    %shading interp;
    % This sets the limits of the colorbar to manual for the first plot
    %caxis manual
    %caxis([-1 1]);
    title([lambdas{p}, ' = ', num2str(lambda(p))],'fontsize',12)
    %title(['\lambda =' num2str(lambda(p))])
    axis off
    
    subplot(2,m+1,m+1+p);
    u = U(:,N-p+1);
    if(u(1)^2>0)
        u = u.*sign(u(1));
    end
    stem(u,'b','linewidth',1);
    %shading interp;
    % This sets the limits of the colorbar to manual for the first plot
    %caxis manual
    %caxis([-1 1]);
    
    title([lambdas{N-p+1}, ' = ', num2str(lambda(N-p+1))],'fontsize',12)
    axis off
    
    
end
subplot(2,m+1,m+1);
    u = U(:,m+1);
    if(u(1)^2>0)
        u = u.*sign(u(1));
    end
    stem(u,'b','linewidth',1);
    %shading interp;
    % This sets the limits of the colorbar to manual for the first plot
    %caxis manual
    %caxis([-1 1]);
    %tt = sprintf('\lambda_{%d}',p);
    title([lambdas{m+1}, ' = ', num2str(lambda(m+1))],'fontsize',12)
    axis off
set(gcf, 'Position',  [200, 200, 1000, 200])
end

function out = lambda_title(N)
out = cell(N,1);
out{1} = '\lambda_1';
out{2} = '\lambda_2';
out{3} = '\lambda_3';
out{4} = '\lambda_4';
out{5} = '\lambda_5';
out{6} = '\lambda_6';
out{7} = '\lambda_7';
out{8} = '\lambda_8';
out{9} = '\lambda_9';
out{10} = '\lambda_{10}';
out{11} = '\lambda_{11}';

end



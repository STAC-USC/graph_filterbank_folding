% Author - Eduardo Pavez <eduardo.pavez.carvelli@gmail.com, pavezcar@usc.edu>
% Copyright Eduardo Pavez,  University of Southern California
%Los Angeles, USA, 02/11/2023
% E. Pavez, B. Girault, A. Ortega, and P. A. Chou. 
%"Two channel filter banks on arbitrary graphs with positive semi-definite variation operators". 
%IEEE Transactions on Signal Processing 2023
%https://arxiv.org/abs/2203.02858
%Example in figure 7
clear;
%image size
m = 8;
n = m;
N = m*n;
%create image point cloud

V = zeros(m*n,2);
k = 1;
for j=1:n
    for i=1:m
        V(k,:) = [i j];
        k = k+1;
    end
end
[E4, D4] = get_image_graph(V,4);%4 connected grid
[E8, D8] = get_image_graph(V,8);% 8 connected grid

W4 = sparse(E4(:,1), E4(:,2), 1, N,N);%sparse adjacency matrix
W8 = sparse(E8(:,1), E8(:,2), 1, N,N);
Wq = W8 - W4;%adjacency matrix with diagonal connections only

[Lists4,partition_4] = bipartition_4_conn(m,n);%vertex partition



%
S = Lists4{1};
R = Lists4{2};
D = diag(sum(W8));
L = D - W8;
Q = L;
Q(S,R) = 0;
Q(R,S) = 0;

B = -L;
B(R,R)=0;
B(S,S)=0;

Btilde = Q\B;
Z = speye(n*m) - Btilde;

Lrw = D\L;

Wbi = W4;
Dbi = diag(sum(Wbi));
Dbi_inv = diag(1./sum(Wbi));

Pbi = Dbi_inv*Wbi;

Pbd = Q\Dbi;

%% plot for bipartite  part
%colormap that makes white=0
cmap = colormap('bone');
cmap = flipud(cmap);
nA = length(S);
nB = length(R);
% relabel verties so that first nA correspond to set S and rest
% correspond to set R
Pbi2 = zeros(N,N);
Pbi2(1:nA,nA+1:end)= Pbi(S,R);
Pbi2(nA+1:end,1:nA)= Pbi(R,S);

%display matrix
imagesc(Pbi2);
colormap(cmap);

%draw lines to mark  vertex partition
line([(N+1)/2, (N+1)/2], [0, N+1],'color','k','linestyle','--');
line([0, N+1],[(N+1)/2, (N+1)/2],'color','k','linestyle','--');

axis square;
caxis([0 1]);
cc=colorbar('location','southoutside');
set(gca,'xtick',[]);
set(gca,'ytick',[])
fnameout = ['fundamental_matrix_image_Pbiv2.png'];
cc.FontSize = 20;
%saveas(gcf,fnameout,'png');

%% plot for block diagonal part
figure;

Pbd2 = zeros(N,N);
Pbd2(1:nA,1:nA)= Pbd(S,S);
Pbd2(nA+1:end,nA+1:end)= Pbd(R,R);
imagesc(Pbd2);
line([(N+1)/2, (N+1)/2], [0, N+1],'color','k','linestyle','--');
line([0, N+1],[(N+1)/2, (N+1)/2],'color','k','linestyle','--');
colormap(cmap);
caxis([0 1]);

axis square;
colormap(cmap);

cc=colorbar('location','southoutside');
set(gca,'xtick',[]);
set(gca,'ytick',[])
fnameout = ['fundamental_matrix_image_Pbdv2.png'];
cc.FontSize = 20;
%saveas(gcf,fnameout,'png');
%% plot for sparsity of fundamental matrix (I- Z)
figure;
PP2 = Pbd2*Pbi2;


imagesc(PP2);
line([(N+1)/2, (N+1)/2], [0, N+1],'color','k','linestyle','--');
line([0, N+1],[(N+1)/2, (N+1)/2],'color','k','linestyle','--');
caxis([0 1]);
axis square;
colormap(cmap);

cc=colorbar('location','southoutside');
set(gca,'xtick',[]);
set(gca,'ytick',[])
fnameout = ['fundamental_matrix_image_PPv2.png'];
cc.FontSize = 20;
%saveas(gcf,fnameout,'png');
%% plot for original adjacency matrix (W)
figure;
W8v2 = zeros(N,N);
W8v2(1:nA,1:nA)= W8(S,S);
W8v2(nA+1:end,nA+1:end)= W8(R,R);

W8v2(1:nA,nA+1:end)= W8(S,R);
W8v2(nA+1:end,1:nA)= W8(R,S);

imagesc(W8v2);
line([(N+1)/2, (N+1)/2], [0, N+1],'color','k','linestyle','--');
line([0, N+1],[(N+1)/2, (N+1)/2],'color','k','linestyle','--');
caxis([0 1]);
axis square;
colormap(cmap);

cc=colorbar('location','southoutside');
set(gca,'xtick',[]);
set(gca,'ytick',[])
fnameout = ['fundamental_matrix_image_W8v2.png'];
cc.FontSize = 20;
%saveas(gcf,fnameout,'png');



function [Lists,partition] = bipartition_vertical(m,n)

partition = zeros(m,n);

for i=1:n
    if(mod(i,2))
        partition(:,i) = ones(m,1);
    end
end
flag = partition(:);

Lists{1} = find(flag==1);
Lists{2} = find(flag==0);
end
function [Lists,partition] = bipartition_horizontal(m,n)
partition = zeros(m,n);

for i=1:m
    if(mod(i,2))
        partition(i,:) = ones(1,n);
    end
end
flag = partition(:);

Lists{1} = find(flag==1);
Lists{2} = find(flag==0);

end
function [Lists4,partition_4] = bipartition_4_conn(m,n)
partition_4 = zeros(m,n);
for i=1:m
    for j=1:n
        if (mod(i+j,2))
            partition_4(i,j)=1;
        end
    end
end
flag = partition_4(:);

Lists4{1} = find(flag==1);
Lists4{2} = find(flag==0);

end

function [edge, distance] = get_image_graph(V,type)
tol = 0.00001;

switch type
    case 4
        rad = 1 + tol;
    case 8
        rad = sqrt(2) + tol;
end

[idx4, d4] = rangesearch(V,V,rad);

E4 = [];
D4 = [];

for i=1:length(idx4)
   ni = idx4{i}(2:end);
   di = d4{i}(2:end);
   sni = length(ni);
   E4 = [E4 ; [repmat(i,sni,1), ni' ]];
   D4 = [D4 ; di' ];
end
E4 = [E4 ; [E4(:,2),E4(:,1)]];
D4 = [D4 ; D4];
[edge, ia, ~] = unique(E4,'rows');
distance = D4(ia);

end

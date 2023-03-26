% Authors - Eduardo Pavez <eduardo.pavez.carvelli@gmail.com, pavezcar@usc.edu>
% Copyright Eduardo Pavez,  University of Southern California
%Los Angeles, USA, 02/11/2023
% E. Pavez, B. Girault, A. Ortega, and P. A. Chou.
%"Two channel filter banks on arbitrary graphs with positive semi-definite variation operators".
%IEEE Transactions on Signal Processing 2023
%https://arxiv.org/abs/2203.02858


%minesotta graph example, requires GRASP toolbox https://github.com/GraSP-toolbox/GraSP
%compare random partitioning with proposed max-cut partitioning

clear;
g = grasp_minnesota();
L = grasp_laplacian_standard(g);



%% this computes the histogram of condition numbers
rng(0);
[condQrandom,condQopt,nonzQrandom,nonzQopt] =compare_condition_number(L);


%% this visualizes an instance of random partition, and the corresponding Q
%matrices 
[A,B,Aran,Bran] = visualize_partition_Q(g);

function [A,B,Aran,Bran] = visualize_partition_Q(g)
L = grasp_laplacian_standard(g);
N = size(L,1);

%compute  maxcut partition
[setA,setB] = spectral_graph_partitioning(L,1);
A = find(setA);
B = find(setB);

%compute random partition
setAran = (rand(N,1)>0.5);

if(nnz(setAran)<N/2)
   setBran = setAran;
   setAran = ~setBran;
else
   setBran = ~setAran; 
end
%
Aran = find(setAran);
Bran = find(setBran);
%display minnesota graph
nodesize = 20;
edgethick = 1.5;
figure;
grasp_show_graph(gca,g,'node_display_size',nodesize,'edge_thickness',edgethick,'background',[]);
title('Minessota Graph');
%saveas(gcf,'minnesota.png','png')
%display partition on minnesota graph

figure;
grasp_show_graph(gca,g,'node_values',(setAran-setBran),'node_display_size',nodesize,'edge_thickness',edgethick,'background',[]);
title('Random Partition on Minessota Graph');
%saveas(gcf,'minnesota_rand.png','png')
figure;
grasp_show_graph(gca,g,'node_values',(setA-setB),'node_display_size',nodesize,'edge_thickness',edgethick,'background',[]);
title('Proposed Max-Cut Partition on Minessota Graph');
%saveas(gcf,'minnesota_maxcut.png','png')
%display same partitions but now on graph induced by Q
Q=L;
Q(Aran,Bran)=0;
Q(Bran,Aran)=0;
g.A_layout = spones( diag(diag(Q)) - Q);
figure;
grasp_show_graph(gca,g,'node_values',(setAran-setBran),'node_display_size',nodesize,'edge_thickness',edgethick,'background',[]);
title('Random Partition on graph of W^{bd}');
%saveas(gcf,'minnesota_rand_Q.png','png')

Q=L;
Q(A,B)=0;
Q(B,A)=0;
g.A_layout = spones( diag(diag(Q)) - Q);
figure;
grasp_show_graph(gca,g,'node_values',(setA-setB),'node_display_size',nodesize,'edge_thickness',edgethick,'background',[]);
title('Proposed Max-Cut Partition on graph of W^{bd}');
%saveas(gcf,'minnesota_msxcut_Q.png','png')
end
function [condQrandom,condQopt,nonzQrandom,nonzQopt] = compare_condition_number(L)
%function to compute condition numbers of various Q matrices generated from
%random sampling and maxcut partitioning
nT = 1000; %numer of realizations
N = size(L,1);
condQrandom = zeros(nT,1);
deg = diag(L);
condD = max(deg)/min(deg);
nonzQrandom = condQrandom;
for i=1:nT
    setAran = (rand(N,1)>0.5);
    if(nnz(setAran)<N/2)
        setBran = setAran;
        setAran = ~setBran;
    else
        setBran = ~setAran;
    end
    A = find(setAran);
    B = find(setBran);
    Q = L;
    Q(A, B) = 0;
    Q(B, A) = 0;
    condQrandom(i) = condest(Q)/condD;
    nonzQrandom(i) = (nnz(Q\L) - N)/(nnz(L)-N);
end

%now compute condition number for maxcut partitioning
[setA,setB] = spectral_graph_partitioning(L,1);
A = find(setA);
B = find(setB);
Q = L;
Q(A, B) = 0;
Q(B, A) = 0;

condQopt = condest(Q)/condD;
nonzQopt = (nnz(Q\L) - N)/(nnz(L)-N);
h=histogram(condQrandom,'Normalization','probability');
hold on;
xline(condQopt,'linewidth',1.5,'color','red');
%histogram(condQopt,'Normalization','probability');
morebins(h);
morebins(h);
xlabel('\kappa(Q)/\kappa(D)','fontsize',16);
ylabel('Probability','fontsize',16);
legend('random partition', 'max-cut partition','fontsize',16,'location','northeast');
%saveas(gcf, 'minnesota_randQ_vs_maxcutQ.png','png');
end
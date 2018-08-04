%% 介绍
%   这是最基本的模型，随机生成点的坐标后，计算其距离矩阵，加上噪声后重构，观察
% 是否与原图形重合。

%% 代码
N = 300;
p = 3;
noise = 0.01;
x = randn(p,N);
x = x - repmat(mean(x,2),[1,N]);

D = pdist2(x',x');
D = D + randn(N,N)*noise;

Dij = D.*D;
Dxj = repmat(sum(Dij,1)./N,[N,1]);
Dix = repmat(sum(Dij,2)./N,[1,N]);
Dxx = ones(N)*sum(Dij(:))/N^2;
B = (Dxj+Dix-Dxx-Dij).*0.5;
[U,S,V]=svd(B);
rex = sqrt(S(1:p,1:p))*U(:,1:p)';
reD = pdist2(rex',rex');

%% 结论
fprintf('the diff between new D and D is %f\n',max(reD(:)-D(:)));
s = sort(diag(S),'descend');
fprintf('the coef of S is %f\n',s(1:p+2));



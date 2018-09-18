function [ubar,kbar] = optimukfun()

% rng('default');
A  = [-1,0;0,-1];
b  = [0;0];
x0 = [50;5];
[X,Y] = ndgrid(25:5:100,0.1:1:5);
W = [X(:),Y(:)];
tpoints = CustomStartPointSet(W);
options = optimoptions('fmincon','Algorithm','interior-point');
problem = createOptimProblem('fmincon','objective',@myfun,'Aineq',A,'bineq',b,'x0',x0,'nonlcon',@mycon,'options',options);

% Construct a GlobalSearch object
% gs = GlobalSearch('Display','iter','StartPointsToRun','bounds-ineqs');
% Construct a MultiStart object based on our GlobalSearch attributes
ms = MultiStart('Display','iter','StartPointsToRun','bounds-ineqs');

% Run GlobalSearch
tic
[x1,myval1] = run(ms,problem,tpoints);
toc

% % Run MultiStart with 15 randomly generated points
% tic
% [x2,myval2] = run(ms,problem,15);
% toc

ubar = x1(1,1);
kbar = x1(2,1);
% ubar = x2(1,1);
% kbar = x2(2,1);

x1
myval1
% x2
% myval2

end 
clear all;
close all;

% create the mesh
Nelems=11;
[X,Y]=meshgrid(0:1/Nelems:1);
T=delaunay(X(:),Y(:));
p=[X(:)';Y(:)'];
t=T';

% assemble stiffness matrix and mass matrix
A=bilin_assembly(@(u,v,du,dv,h)du{1}.*dv{1}+du{2}.*dv{2},p,t);
M=bilin_assembly(@(u,v,du,dv,h)u.*v,p,t);

% find boundary indices
D=find((p(1,:)==0.0)+(p(1,:)==1.0)+(p(2,:)==1.0)+(p(2,:)==0.0));
I=setdiff(1:size(p,2),D);

% initialize and compute the solution of the equation -\delta u = f
% with unit load f=1
u=zeros(size(p,2),1);
f=ones(size(p,2),1);
u(I)=A(I,I)\(M(I,I)*f(I));

% visualize solution
figure;
trisurf(t', X(:), Y(:), u);

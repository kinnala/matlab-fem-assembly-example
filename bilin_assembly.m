% Simplified P1 assembly. Assembles a given bilinear form
% 
% CALLING SYNTAX IS 
%
% function B=bilin_assembly(bilin,p,t)
%
%  bilin = function handle to integrand of the bilinear form
%          given as bilin(u,v,du,dv,h), where
%          u  - values of function u
%          v  - values of function v
%          du{1} - x-derivative of function u
%          du{2} - y-derivative of function u
%          dv{1} - x-derivative of function v
%          dv{2} - y-derivative of function v
%          h - the mesh parameter
%  p, t  = the mesh on which the matrix will be assembled,
%          each column of p is one vertex and each column
%          of t is one element, three vertex indices per column.
%
%  B     = matrix related to bilinear form bilin

function B=bilin_assembly(bilin,p,t)
    nv=size(p,2);
    nt=size(t,2);
    
    B=sparse(nv,nv);

    % quadrature rule
    qp=[1/6,2/3,1/6;2/3,1/6,1/6];
    qw=[1/6,1/6,1/6]';
    
    % basis functions at quadrature points
    phi{1}=1-qp(1,:)-qp(2,:);
    phi{2}=qp(1,:);
    phi{3}=qp(2,:);

    % basis function gradients at quadrature points
    gradhat_phi{1}=repmat([-1;-1],1,length(qw));
    gradhat_phi{2}=repmat([1;0],1,length(qw));
    gradhat_phi{3}=repmat([0;1],1,length(qw));

    % affine mappings to all triangles
    A=cell(2,2);
    A{1,1}=p(1,t(2,:))'-p(1,t(1,:))';
    A{1,2}=p(1,t(3,:))'-p(1,t(1,:))';
    A{2,1}=p(2,t(2,:))'-p(2,t(1,:))';
    A{2,2}=p(2,t(3,:))'-p(2,t(1,:))';

    % determinants of all affine mappings
    detA=A{1,1}.*A{2,2}-A{1,2}.*A{2,1};

    % all inverse mappings
    invAt=cell(2,2);
    invAt{1,1}=A{2,2}./detA;
    invAt{2,1}=-A{1,2}./detA;
    invAt{1,2}=-A{2,1}./detA;
    invAt{2,2}=A{1,1}./detA;

    % find mesh parameter
    h = sqrt(abs(detA))

    for j=1:3
        % function values
        u=repmat(phi{j},nt,1);
        du={};
        du{1}=bsxfun(@times,invAt{1,1},gradhat_phi{j}(1,:))+bsxfun(@times,invAt{1,2},gradhat_phi{j}(2,:));
        du{2}=bsxfun(@times,invAt{2,1},gradhat_phi{j}(1,:))+bsxfun(@times,invAt{2,2},gradhat_phi{j}(2,:));
        
        for i=1:3
            v=repmat(phi{i},nt,1);
            dv={};
            dv{1}=bsxfun(@times,invAt{1,1},gradhat_phi{i}(1,:))+bsxfun(@times,invAt{1,2},gradhat_phi{i}(2,:));
            dv{2}=bsxfun(@times,invAt{2,1},gradhat_phi{i}(1,:))+bsxfun(@times,invAt{2,2},gradhat_phi{i}(2,:));

            B=B+sparse(t(i,:),t(j,:),(bilin(u,v,du,dv,h)*qw).*abs(detA),nv,nv);
        end
    end
end

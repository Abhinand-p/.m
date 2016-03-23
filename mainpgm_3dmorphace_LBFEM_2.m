
clc;
clear all;
close all;
getd = @(p)path(p,path); 
getd('toolbox_signal/');
getd('toolbox_general/');
getd('toolbox_graph/');
getd('triangulation_toolbox/');
getd('matlabmesh/');

lamdall=[];
lamda10=[];
scoreall = [];
for t=13

cd('/Users/ambikasatish/Desktop/PhD_26.3.15/PhD/database/3D morphace database/Morphace scans');
% cd('sub2_session1');
    str=strcat('10_',int2str(t),'.png');   
    eval('I=imread(str);');

% I = imread('morphace4.png');

cd('/Users/ambikasatish/Desktop/PhD_26.3.15/PhD/MATLAB/Biometrics/LB operator/LB_FEM');

try
    I = rgb2gray(I);
end
 
% I = imcrop(I,[218 284 188 136]);%for t=1
I = imcrop(I,[188 290 176 124]);%for t=13
% I = imcrop(I,[208 286 188 136]);%for t=2
I =imresize(I,[128 128]);
[m n] = size(I);

%%%%%%%%%%%%%%%%% Triangulation %%%%%%%%%%%%%%%%%%%%%

[dt z maxerr] = aIDT(I,0,500);
A = aReCon(dt,z);
figure;imshow(I,[]);xlabel('Original Image');
figure;imshow(A,[]);ylabel('Approximated Image');xlabel(['Absolute Error: ',num2str(maxerr)]);
figure;imshow(A,[]);xlabel('Approximated Image and Delaunay Triangulation');
hold on;triplot(dt,'r');xlabel('Delaunay Triangulation');axis([0 n 0 m]);

%%%%%%%%%%%%%%%%%%%%%%%%%%compute_faces.m file %%%%%%%%%%%%%%%%%%%%%%
[r c] = size(dt);
for i = 1:r
    for j = 1:c
        dt1(i,j) = dt(i,j);
    end 
end


faces = compute_delaunay(dt.X);
faces = faces';

peri.faces = faces;
peri.vertices = dt.X;
peri;
%%%%%%%%%%%%%%%%%%%%%%%FEM [A, C] =FEM(peri);
tri   = peri.faces;
coord = peri.vertices;
n = size(coord,1);
%%%%%%%Adjacency matrix
%%%%%%% The incidence matrix in the old code is the adjaency matrix

Adj = sparse(tri(:,[1 2 3]),tri(:,[2 3 1]),1,n,n);
Adj = double(Adj|Adj'); % adjacent nodes
[i,j] = find(Adj);  
[ei,ej] = find(Adj(:,i) & Adj(:,j)); % common adjacent triangle-edge
 e1 = find([1; diff(ej)]); % 1:2:end  % ->ignors singularities
e2 = e1+1;                % 2:2:end

%%%%%%%%%%%% ------------------------------
%%%%%%%%%%%% coordiantes of adjacent nodes and of common adjacent triangle-edge
pi = coord(i,:);
pj = coord(j,:);
qi = coord(ei(e1),:); % qi = coord(ei(1:2:end),:);
qj = coord(ei(e2-1),:); 
% qj = coord(ei(2:2:end),:);

%%%%%%%%% distances
dii = pi-qi; dij = pi-qj;
dji = pj-qi; djj = pj-qj;

norm = @(x) sqrt(sum(x.^2,2)); % norm of a vector
cross2 = @(x,y) x(1)*y(2) - x(2)*y(1);

%%%%%%%% ------------------------------
%%%%%%%% A matrix - area matrix
area = @(vi,vj) norm(cross2(vi,vj))./2;
A = sparse(i,j,(area(dii,dji) + area(dij,djj))./12,n,n);
A = A + diag(sum(A,2));


% % % % % % % ------------------------------
% % % % % % % C matrix - cotangent matrix
cotan = @(vi,vj) cot( acos(dot(vi,vj,2)./(norm(vi).*norm(vj))) );
C = sparse(i,j,-(cotan(dii,dji) + cotan(dij,djj))./2,n,n);
C = C - diag(sum(C,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[V, D] = eigs(C,A,50,'sm');
value = sort(diag(D),'descend');
lamda = value(1:10,1);
lamdall = [lamdall,value];
lamda10 = [lamda10,lamda];
lamdaexp = lamda*10000;
lamdalog = log10(lamdaexp);
score = mean(lamdalog);
scoreall = [scoreall,score];
end


% % % % % % % % % % % % % %%%%%%%%% TO compute GPS %%%%%%%%%%%%% to be corrected
% % % % % % % % % % % % % GPS41 = [];
% % % % % % % % % % % % % for i = 1:size(D,1)
% % % % % % % % % % % % %     phi = V(:,i);
% % % % % % % % % % % % %     la = diag(D);
% % % % % % % % % % % % %     GPS1 = phi./sqrt(la(i));
% % % % % % % % % % % % %     GPS41 = [GPS41, GPS1];
% % % % % % % % % % % % % end
% % % % % % % imshow(GPS2);
% end
% % peri1 = lb_smooth([],peri,0.2,100,V,D);
% % figure_patch(peri1,[0.7 0.7 0.6],0.5)
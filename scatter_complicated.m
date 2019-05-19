% March, 2012
%
% This program sets up a simple external scattering problem for Helmholtz.


function main_develop

clear;
close all;
clc;

addpath('~/Box Sync/research/code/export_fig/');

grid.name = 'contour and disc geometry structure';
tree.name = 'tree structure for storing ';
param = LOCAL_get_parameters();
grid = LOCAL_get_scatterer_geometry(grid,param);
% LOCAL_plot_geometry(grid, param, false);
tree = LOCAL_generate_index_tree(tree,param);
tree = LOCAL_generate_neighbour_list(tree,grid,param);
tic;
tree = LOCAL_compress_all_levels(tree,grid,param);
tree = LOCAL_generate_scattering_matrices(tree);
fprintf('Time to compute scattering matrix for complicated geometery = %0.4f\n', toc);
tree = LOCAL_invert_HBS_matrix(tree);
return


function tree = LOCAL_invert_HBS_matrix(tree)
    for nid=tree.nnodes:-1:2
        if tree.node{nid}.isleaf
            X = tree.node{nid}.A;
        else
            lid = 2*nid;
            rid = 2*nid+1;
            At_lr = tree.node{nid}.At_lr;
            At_rl = tree.node{nid}.At_rl;            
            X1 = tree.node{lid}.D;
            X2 = tree.node{rid}.D;
            X = [X1, At_lr; At_rl, X2];
        end
        
        Xinv = inv(X);
        T = tree.node{nid}.T;
        k = tree.node{nid}.k;
        Ustar = [eye(k), T];
        D = inv(Ustar*Xinv*Ustar');
        Xinv_U = Xinv*Ustar';
        Ustar_Xinv = Ustar*Xinv;
        E = Xinv_U*D;
        Fstar = D*Ustar_Xinv;
        tree.node{nid}.G = X - E*Ustar_Xinv;
        tree.node{nid}.E = E;
        tree.node{nid}.Fstar = Fstar;
        tree.node{nid}.D = D;
    end
    
    G = [tree.node{2}.D, tree.node{1}.At_lr; tree.node{1}.At_rl, tree.node{3}.D];
    tree.node{1}.G = inv(G);
return

function tree = LOCAL_compress_all_levels(tree,grid,param)
    for nid=tree.nnodes:-1:1
        if tree.node{nid}.isleaf
            indskel = tree.node{nid}.itau;
        else
            % concatenate indskels of left and right children
            lid = tree.node{nid}.left_child;
            rid = tree.node{nid}.right_child;
            indskel = [tree.node{lid}.Ik, tree.node{rid}.Ik];
        end
        
        % get the geometry at current node
        C = grid.C(:, sort(indskel));
        % get the circumcircle
        [xc,yc,R,~] = LOCAL_get_circumcircle(C);
        % get proxy circle
        Cproxy = LOCAL_get_proxy_circle(xc, yc, param.proxyscale*R, grid, param);
        % loop over the neighbour list to find points inside the proxy
        % circle
        inside = [];        
        for neid=1:length(tree.node{nid}.neighbours)
            nei = tree.node{nid}.neighbours(neid);
            if tree.node{nei}.isleaf
                inside = [inside, tree.node{nei}.itau];
            else
                lneid = tree.node{nei}.left_child;
                rneid = tree.node{nei}.right_child;
                inside = [inside, tree.node{lneid}.Ik, tree.node{rneid}.Ik];
            end            
        end
        % find which indices are inside the proxy circle
        inside = inside( (grid.C(1,inside) - xc).^2  + (grid.C(4,inside)-yc).^2 < (param.proxyscale*R)^2);
        
%         
        Aout = [LOCAL_get_A_single_offd(grid.C,  inside,  indskel, param);...
                    LOCAL_get_A_sing(Cproxy, 1:param.M, grid.C, indskel, param)];
        Ain = [LOCAL_get_A_single_offd(grid.C,  indskel, inside, param),...
                    LOCAL_get_A_sing(grid.C, indskel, Cproxy, 1:param.M, param)];

%         Aout = [LOCAL_get_A_single_offd(grid.C,  inside,  indskel, param);...
%                     LOCAL_outgoing_to_proxy(Cproxy, 1:param.M, grid.C, indskel, param)];
%         Ain = [LOCAL_get_A_single_offd(grid.C,  indskel, inside, param),...
%                     LOCAL_incoming_from_proxy(grid.C, indskel, Cproxy, 1:param.M, param)];
%                                     
        % Compute the skeletons.
       [T,I] = id_decomp([Aout;Ain'], param.acc);
       % get rank
       k  = size(T,1);
       
       tree.node{nid}.k = k;
       tree.node{nid}.T = T;
       tree.node{nid}.Ik = indskel(I(1:k));
       tree.node{nid}.I = I;       
    end
    
    
    % construct self interation matrices    
    for nid=tree.nnodes:-1:1
        if tree.node{nid}.isleaf
            itau = tree.node{nid}.itau;
            tree.node{nid}.A = LOCAL_get_A_single_diag(grid.C,itau,param);                            
        end
    end
            
    
    % construct sibling interation matrices
    for nid=tree.nnodes:-1:1
        if ~tree.node{nid}.isleaf
            lid = 2*nid;
            rid = 2*nid+1;
            tree.node{nid}.At_lr = LOCAL_get_A_single_offd(grid.C, tree.node{lid}.Ik , tree.node{rid}.Ik, param);
            tree.node{nid}.At_rl = LOCAL_get_A_single_offd(grid.C, tree.node{rid}.Ik , tree.node{lid}.Ik, param);
        end
    end
            
return

function [xc,yc,R,Z] = LOCAL_get_circumcircle(C)
    xc = mean(C(1,:));
    yc = mean(C(4,:));
    dist = sqrt((C(1,:)-xc).^2 + (C(4,:)-yc).^2);
    R = max(dist(:));    
    nt = 100;
    h = 2*pi/nt;
    theta = 0:h:2*pi-h;
    Z = R*exp(sqrt(-1)*theta);
    Z = [real(Z);imag(Z)] + [xc;yc];
return

function tree = LOCAL_generate_neighbour_list(tree,grid,param)
%     tree.node{1}.neighbours = [];
    for nid=1:tree.nnodes
        % get relevant part of the geometry 
        C = grid.C(:,tree.node{nid}.itau);
        % get circumcircle of the geometry
        [xc,yc,R,~] = LOCAL_get_circumcircle(C);                        
        Rproxy = R*param.proxyscale;
        % find out if any neighbours on the current level have points
        % inside the proxy circle
        tree.node{nid}.neighbours = [];
        for jid=1:tree.nnodes
            % check if on the same level as nid
            if tree.node{jid}.level == tree.node{nid}.level && nid~=jid
                % check if any points of jid inside nid's proxy surface
                Cjid = grid.C(:, tree.node{jid}.itau);
                dist = (Cjid(1,:)-xc).^2 + (Cjid(4,:)-yc).^2;
                if any(dist < Rproxy^2)
                    tree.node{nid}.neighbours(end+1) = jid;
                end
            end
        end
        % plot them       
%         scatter(C(1,:), C(4,:)); hold on; scatter(Z1(1,:), Z1(2,:));
%         scatter(Z2(1,:), Z2(2,:)); hold on; scatter(grid.C(1, tree.node{3}.itau), grid.C(4,tree.node{3}.itau));
    end
        
return

function tree = LOCAL_generate_scattering_matrices(tree)
    for nid=tree.nnodes:-1:1
        if tree.node{nid}.isleaf
            k = tree.node{nid}.k;
            T = tree.node{nid}.T;
            Ustar = [eye(k), T];
            tree.node{nid}.S = Ustar * inv(tree.node{nid}.A) * Ustar';
        else
            tree = LOCAL_form_scattering_matrix_non_leaf_node(nid, tree);
        end
    end
return

function tree = LOCAL_form_scattering_matrix_non_leaf_node(nid, tree)
    k = tree.node{nid}.k;
    T = tree.node{nid}.T;
    Ustar = [eye(k), T];
    lid = 2*nid;
    rid = 2*nid+1;
    Sl = tree.node{lid}.S;
    Sr = tree.node{rid}.S;
    At_lr = tree.node{nid}.At_lr;
    At_rl = tree.node{nid}.At_rl;
    Z = zeros(size(Sl,1), size(Sr,2));
    tree.node{nid}.S = Ustar * inv([eye(size(Sl)), Sl*At_lr;  ...
                                 Sr*At_rl, eye(size(Sr))]) ...
                               *[Sl, Z;     ...
                                 Z', Sr] * Ustar';
return

function tree = LOCAL_generate_index_tree(tree,param)
    if param.ng ~= 2^(nextpow2(param.ng))
        warning('number of contours must be powers of 2, tree may not be created');
        tree = 0;
        return;        
    end
    % root is all geometry
    nlevels = log2(param.ng)+1;
    tree.nnodes = 2^nlevels - 1;
    id = 1:param.ng;
    tree.node{1}.gid = id;    
    if tree.nnodes > 1
        tree.node{1}.left_child = 2;
        tree.node{1}.right_child = 3;
        tree.node{1}.isleaf = false;
    else
        tree.node{1}.left_child = nan;
        tree.node{2}.right_child = nan;
        tree.node{1}.isleaf = true;
    end
    tree.node{1}.level = 0;
    tree.node{1}.parent = nan;
    tree.node{1}.itau = 1:param.ng*param.N;
    tree.nlevels = nlevels;  
    N = param.N;
    for l=1:nlevels-1
        % split the array (2^l) ways 
        nsplit = 2^l;
        k = 1;
        out = mat2cell(id,1,(param.ng/nsplit)*ones(1,nsplit));    
        for i=2^l:2^(l+1)-1
            tree.node{i}.gid = out{k};
            gid = tree.node{i}.gid;
            tree.node{i}.level = l;
            tree.node{i}.parent = floor(i/2);            
            if l==nlevels-1
                tree.node{i}.isleaf = true;
                tree.node{i}.left_child = nan;
                tree.node{i}.right_child = nan;                 
                
                tree.node{i}.indskel = (gid-1)*N + 1 : gid*N;
                tree.node{i}.itau = tree.node{i}.indskel;
            else
                tree.node{i}.isleaf = false;
                tree.node{i}.left_child = 2*i;
                tree.node{i}.right_child = 2*i+1;                
                tree.node{i}.itau = (gid(1)-1)*N + 1 : gid(end)*N;
            end            
            k = k + 1; 
        end        
    end    
return

function param = LOCAL_get_parameters()
    % number of stars in each direction
    param.ngx = 2;
    param.ngy = 2;
    param.ng = param.ngx * param.ngy;
    % number of grid points in a star
    param.N = 200;
    param.h = 2*pi/param.N;
    % number of points on the circular proxy surface
    param.M = 50;
    % spacing of the stars
    param.spx = 4;
    param.spy = 3;
    % wave number in helmholtz operator
    param.kh = 2;
    % rank of off-diagonal blocks
    param.k = param.N/2;
    param.acc = 1e-10;
    % proxy circle scale
    param.proxyscale = 1.5;
return

function grid = LOCAL_get_scatterer_geometry(grid,param)
    grid.N = param.N;
    grid.M = param.M;
    ngx = param.ngx;
    ngy = param.ngy;
    ng = param.ng;
    spx = param.spx;
    spy = param.spy;
    N = param.N;
    % create an array of centers for the stars            
    [grid.xc, grid.yc] = meshgrid(0:spx:ngx*spx-1, 0:spy:ngy*spy-1);
    grid.xc = grid.xc + linspace(0,12,size(grid.xc,1))';
    grid.xc = grid.xc(:);
    grid.yc = grid.yc(:);
    % matrix to store coordinates of geomtries 
    grid.C = zeros(6, N*ng);
    fprintf('generating star geometry grid\n');
    for i=1:ngx
        for j=1:ngy
            k = (j-1 + (i-1)*ngy);            
            grid.C(:, N*k+1:N*(k+1)) = LOCAL_get_geometry_star(N, [grid.xc(k+1); grid.yc(k+1)]);            
        end
    end
    grid.diameter = 2*max(sqrt(grid.C(1,1:N).^2 + grid.C(4,1:N).^2));  % The diameter of the contour.
return

function Z = LOCAL_get_proxy_circle(xc, yc, R, grid, param)
    M = param.M;
    h = 2*pi/M;
    tt = 0:h:2*pi-h;
    Z  =     [xc + R*cos(tt);...
                 - R*sin(tt);...
                 - R*cos(tt);...
              yc + R*sin(tt);...
                 + R*cos(tt);...
                 - R*sin(tt)];
return

function LOCAL_plot_geometry(grid,param,save)
    N = grid.N;
    figure;
    for id=1:param.ng
        k = id-1;
        Cid = grid.C(:,N*k+1:N*(k+1));
        plot( Cid(1,:), Cid(4,:), '-k','Linewidth',2); hold on;               
    end
    hold off;
    axis tight;
    axis equal;
    set(gca,'XColor','none'); set(gca,'YColor','none')
    if save
        export_fig('figures/geometry1.pdf', '-q101', '-pdf', '-transparent');
    end
return

function A = LOCAL_get_A_sing(C1,ind1,C2,ind2,param)
    kh = param.kh;
    h  = param.h;

    [Y_g1,  X_g1  ] = meshgrid(C2(1,ind2), C1(1,ind1));
    [Y_g2,  X_g2  ] = meshgrid(C2(4,ind2), C1(4,ind1));
    [Y_dg1, ~] = meshgrid(C2(2,ind2), C1(2,ind1));
    [Y_dg2, ~] = meshgrid(C2(5,ind2), C1(5,ind1));    
    ima = sqrt(-1);
    nn1 = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    nn2 = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
    dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
    A   = ((nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd)) + ...
          ima*kh*besselh(0, kh*dd)).*h.*sqrt(Y_dg1.^2 + Y_dg2.^2);    
return
       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes diagonal blocks of the coefficient matrix.
% Note that this should be used ONLY for an index set ind that is
% entirely contained within one contour.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = LOCAL_get_A_single_diag(C,ind,param)

kh       = param.kh;
ntot     = size(C,2);
h        = param.h;
nloc     = length(ind);

[Y_g1,   X_g1   ] = meshgrid(C(1,ind), C(1,ind));
[Y_g2,   X_g2   ] = meshgrid(C(4,ind), C(4,ind));
[Y_dg1,  ~  ] = meshgrid(C(2,ind), C(2,ind));
[Y_dg2,  ~  ] = meshgrid(C(5,ind), C(5,ind));

ima = sqrt(-1);
nn1 = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
nn2 = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2) + eye(nloc);
A   = ((nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd)) + ...
       ima*kh*besselh(0, kh*dd)).*h.*sqrt(Y_dg1.^2 + Y_dg2.^2);
MU6 = [ 0.2051970990601250e1 + 0.2915391987686505e1;...
       -0.7407035584542865e1 - 0.8797979464048396e1;...
        0.1219590847580216e2 + 0.1365562914252423e2;...
       -0.1064623987147282e2 - 0.1157975479644601e2;...
        0.4799117710681772e1 + 0.5130987287355766e1;...
       -0.8837770983721025   - 0.9342187797694916];
D = abs(ind'*ones(1,length(ind)) - ones(length(ind),1)*ind);
D = abs(mod(D + 6, ntot) - 6);
[ii1, ii2] = find(D.*(D<7));
for j = 1:length(ii1)
   i1 = ii1(j);
   i2 = ii2(j);
   d  = D(i1,i2);
   A(i1,i2) = (1 + MU6(d))*A(i1,i2);
end
for j = 1:nloc
  A(j,j) = - 2*ima;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the field caused by point charges.
%   xx_trg   = locations of target points.
%   xx_src   = locations of sources points.
%   qq       = the strength of the point charges.
%   kh       = helmholtz parameter.
%
% The output u is the vector
%
%   u(i) = sum_{j} qq(j) * H_0(kh*|xx_trg(:,i) - xx_src(:,j)|)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = LOCAL_pot_from_charges(xx_trg,xx_src,qq,kh)

ntrg  = size(xx_trg,2);
nsrc  = size(xx_src,2);
dist1 = (xx_trg(1,:)' * ones(1,nsrc)) - ones(ntrg,1)*xx_src(1,:);
dist2 = (xx_trg(2,:)' * ones(1,nsrc)) - ones(ntrg,1)*xx_src(2,:);

dist  = sqrt(dist1.*dist1 + dist2.*dist2);

u = besselh(0,kh*dist)*reshape(qq,numel(qq),1);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates a "smooth star" centered at xxc.
% The contour is parameterized as follows:
%
%   C(1,i) =                      x1 coordinate of node i
%   C(2,i) =        derivative of x1 coordinate of node i
%   C(3,i) = second derivative of x1 coordinate of node i
%   C(4,i) =                      x2 coordinate of node i
%   C(5,i) =        derivative of x2 coordinate of node i
%   C(6,i) = second derivative of x2 coordinate of node i
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,curvelen,xx_int,xx_ext] = LOCAL_get_geometry_star(ntot,xxc)

r        = 0.3;
k        = 5;
tt       = linspace(0,2*pi*(1 - 1/ntot),ntot);
C        = zeros(6,ntot );
C(1,:)   =   1.5*cos(tt) + (r/2)*            cos((k+1)*tt) + (r/2)*            cos((k-1)*tt);
C(2,:)   = - 1.5*sin(tt) - (r/2)*(k+1)*      sin((k+1)*tt) - (r/2)*(k-1)*      sin((k-1)*tt);
C(3,:)   = - 1.5*cos(tt) - (r/2)*(k+1)*(k+1)*cos((k+1)*tt) - (r/2)*(k-1)*(k-1)*cos((k-1)*tt);
C(4,:)   =       sin(tt) + (r/2)*            sin((k+1)*tt) - (r/2)*            sin((k-1)*tt);
C(5,:)   =       cos(tt) + (r/2)*(k+1)*      cos((k+1)*tt) - (r/2)*(k-1)*      cos((k-1)*tt);
C(6,:)   = -     sin(tt) - (r/2)*(k+1)*(k+1)*sin((k+1)*tt) + (r/2)*(k-1)*(k-1)*sin((k-1)*tt);
curvelen = 2*pi;

rmin = sqrt(min(C(1,:).^2 + C(4,:).^2));
rmax = sqrt(max(C(1,:).^2 + C(4,:).^2));

%%% Construct interior points.
ttint  = 2*pi*sort(rand(1,3));
xx_int = 0.5*rmin*[cos(ttint); sin(ttint)];

%%% Construct exterior points.
ttext  = 2*pi*sort(rand(1,5));
xx_ext = 1.5*rmax*[cos(ttext); sin(ttext)];

%%% Shift to center at xxc.
C([1,4],:) = C([1,4],:) + xxc*ones(1,size(C,2));
xx_int     = xx_int     + xxc*ones(1,size(xx_int,2));
xx_ext     = xx_ext     + xxc*ones(1,size(xx_ext,2));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute an off-diagonal block of the matrix discretizing a single contour.
% Note that this INCLUDES corrections near the diagonal for handling the
% singularity in the kernel.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = LOCAL_get_A_single_offd(C,ind1,ind2,param)

ntot     = size(C,2);
kh       = param.kh;
h        = param.h;

[Y_g1,  X_g1  ] = meshgrid(C(1,ind2), C(1,ind1));
[Y_g2,  X_g2  ] = meshgrid(C(4,ind2), C(4,ind1));
[Y_dg1, ~ ] = meshgrid(C(2,ind2), C(2,ind1));
[Y_dg2, ~ ] = meshgrid(C(5,ind2), C(5,ind1));

ima = sqrt(-1);
nn1 = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
nn2 = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
A   = ((nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd)) + ...
       ima*kh*besselh(0, kh*dd)).*h.*sqrt(Y_dg1.^2 + Y_dg2.^2);
MU6 = [ 0.2051970990601250e1 + 0.2915391987686505e1;...
       -0.7407035584542865e1 - 0.8797979464048396e1;...
        0.1219590847580216e2 + 0.1365562914252423e2;...
       -0.1064623987147282e2 - 0.1157975479644601e2;...
        0.4799117710681772e1 + 0.5130987287355766e1;...
       -0.8837770983721025   - 0.9342187797694916];
D = abs(ind1'*ones(1,length(ind2)) - ones(length(ind1),1)*ind2);
D = abs(mod(D + 6, ntot) - 6);
[ii1, ii2] = find(D.*(D<7));
for j = 1:length(ii1)
  i1 = ii1(j);
  i2 = ii2(j);
  d  = D(i1,i2);
  A(i1,i2) = (1 + MU6(d))*A(i1,i2);
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes off-diagonal blocks of the coefficient matrix,
% but without adding any quadrature corrections. This is intended to be 
% used for blocks for which ind1 and ind2 belong to different contours.
% (It can also be used for well-separated blocks regardless of location.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A = LOCAL_get_A_offd_noquad(C,ind1,ind2,param)

kh = param.kh;
h  = param.h;

[Y_g1,  X_g1  ] = meshgrid(C(1,ind2), C(1,ind1));
[Y_g2,  X_g2  ] = meshgrid(C(4,ind2), C(4,ind1));
[Y_dg1, ~ ] = meshgrid(C(2,ind2), C(2,ind1));
[Y_dg2, ~ ] = meshgrid(C(5,ind2), C(5,ind1));

ima = sqrt(-1);
nn1 = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
nn2 = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
dd  = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
A   = ((nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd)) + ...
       ima*kh*besselh(0, kh*dd)).*h.*sqrt(Y_dg1.^2 + Y_dg2.^2);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function evaluates a layer potential with the combined field kernel.
% xx     is the matrix of target points, 
% C      is the contour.
% sigma  is the source distribution on C.
% params contains the Helmholtz parameter, the grid spacing "h", etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vv = LOCAL_evalpot(xx, C, sigma, params)

kh = params(2);
h  = params(3);

n = size(C,2);
m = size(xx,2);

X_g1  = xx(1,:)'*ones(1,n);
X_g2  = xx(2,:)'*ones(1,n);
Y_g1  = ones(m,1)*C(1,:);
Y_g2  = ones(m,1)*C(4,:);
Y_dg1 = ones(m,1)*C(2,:);
Y_dg2 = ones(m,1)*C(5,:);

ima  = sqrt(-1);
nn1  = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
nn2  = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
dd   = sqrt((Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2);
EVAL = ((nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2)).*(1./dd).*(-kh*besselh(1, kh*dd)) + ...
        ima*kh*besselh(0, kh*dd)).*h.*sqrt(Y_dg1.^2 + Y_dg2.^2);
   
vv = EVAL*sigma;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the operator which maps charges q on disk to contour.
% X      is the contour.
% Z      is the disc.
% ind1   is the indices to be used in the contour
% ind2   is the indices to be used on the proxy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = LOCAL_incoming_from_proxy(X, ind1, Z, ind2, param)

kh = param.kh;
M = size(Z,2);      
% coordinates of the contour and disk in meshgrid format
[X1,X2] = meshgrid(Z(1,ind2), X(1,ind1));
[Y1,Y2] = meshgrid(Z(2,ind2), X(4,ind1));
% distance matrix
dd = sqrt((X1-X2).^2 + (Y1-Y2).^2);
% C=\phi_k(x_i, z_j) * 2*\pi/M
C = (sqrt(-1)/4) * besselh(0, kh*dd) * 2*pi/M;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the operator which maps charges q on disk to contour.
% X      is the contour.
% Z      is the disc.
% ind1   is the indices to be used in the contour
% ind2   is the indices to be used in the disc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = LOCAL_outgoing_to_proxy(Z, ind2, X, ind1, param)

% get parameters
kh = param.kh;
h = param.h;

% coordinates of the contour and disk in meshgrid format
[X1,X2] = meshgrid(X(1,ind1), Z(1,ind2)); % X1 = x-coordinates of disk, X2 = x-coordinates of contour
[Y1,Y2] = meshgrid(X(4,ind1), Z(4,ind2)); % Y1 = y-coordinates of disk, Y2 = y-coordinates of contour
% distance matrix
dd = sqrt((X1-X2).^2 + (Y1-Y2).^2);

% normals for the contour
[Nx1, ~] = meshgrid(X(2,ind1), Z(2,ind2)); % x-component of normal
[Ny1, ~] = meshgrid(X(5,ind1), Z(5,ind2)); % y-component of normal
abs_nn = sqrt(Nx1.^2 + Ny1.^2);
nn1 = Ny1./abs_nn;
nn2 = -Nx1./abs_nn;

B  = ((nn1.*(X1-X2) + nn2.*(Y1-Y2)).*(1./dd).*(-kh*besselh(1, kh*dd)) + ...
       sqrt(-1)*kh*besselh(0, kh*dd)).*h.*abs_nn;
return

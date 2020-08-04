%PartB load the resource
%%under sampled data
%%
load('SLphan.mat');
im = SLphan;
figure(1); clf
imagesc(im); colormap(gray);title('original image');

%full projection 
angles = [0:179];
g = radon(im,angles);
figure(2);
imagesc(g);colormap gray; title('radon transform' );


%undersampled 
N = 180;
[gx,gy] = size(g);
resizeg = g(1:180,1:180);

%limited angles 
mask = ones(N,N);
mask(1:N,70:110) = 0;
lim = mask.*resizeg;
figure(3);
imagesc(lim); colormap(gray);title('limited angle');

%undersample image 
mask2 = ones(N,N);
for i = 1:2:59
	mask2(1:N,3*i:3*(i+1)) = 0;
end

under = mask2.*resizeg;
figure(4);
imagesc(under); colormap(gray);title('undersampled');

%Using isotropic Laplacian
h = 1/(N+1);
D1 = sparse(N,N);
D1(1,1) = 1;
for k = 2:N
    D1(k-1,k) = -1;
    D1(k,k) = 1;
end

Dx = kron(speye(N),D1);
Dy = kron(D1,speye(N));
Lapl = -(Dx'*Dx + Dy'*Dy);

imx = reshape(Dx*reshape(resizeg,[],1),N,N);
imy = reshape(Dy*reshape(resizeg,[],1),N,N);
imlap = reshape(Lapl*reshape(resizeg,[],1),N,N);

figure(6); clf
subplot(2,2,1);imagesc(resizeg); colormap(gray);title('Original f');
subplot(2,2,2);imagesc(imx); colormap(gray);title('f_x');
subplot(2,2,3);imagesc(imy); colormap(gray);title('f_y');
subplot(2,2,4);imagesc(imlap); colormap(gray);title('\nabla^2 f');


%% now we try PDE with Dirichlet B.c.s
ind2 = find(mask2==1);
outd = setdiff([1:N*N],ind2);
g = zeros(N*N,1);
g(ind2) = resizeg(ind2); % this is an "arbitrary" data vector;


outLapl = Lapl;
outLapl(ind2,:) = 0; % zeros(length(ind2),N*N);

% First try solving the PDE Lap f = g
finP = speye(N*N);
finP(outd,outd) = 0;
f1 = (finP+outLapl)\g;
f1 = reshape(f1,N,N);

figure(8); clf; hold on;
subplot(2,1,1); imagesc(reshape(g,N,N));colormap(gray)
subplot(2,1,2); imagesc(f1);colormap(gray); title('PDE + Dirichlet b.c.s \alpha= ');

%filtered backprojection
gcorrect = zeros(185,180);
gcorrect(1:180,1:180) = f1;
n = 128;
ig0 = iradon(gcorrect,angles,1,n);
figure(7);
subplot(1,2,1);imagesc(im);colormap gray;title('original image' );
subplot(1,2,2);imagesc(ig0);colormap gray;title('filtered back-projection' );


%%
I180 = speye(180*180);
N = 180;
mind = find(mask2==1);
A = I180(mind,:);

g = A*reshape(resizeg,[],1); % this is the reduced size data;
alpha0 = 1;
fra = inv(A'*A + alpha0*I180)*A'*g;
fraim = reshape(fra,N,N);

figure(5);clf;
subplot(2,2,1);imagesc(resizeg); colormap(gray);title('original image');
subplot(2,2,2);imagesc(under); colormap(gray);title('masked image');
subplot(2,2,3);imagesc(fraim); colormap(gray);title(['recon with 0-Tikononv \alpha= ',num2str(alpha0)]);


alphaIso = 1e-3;
Atot = [A;sqrt(alphaIso)*Dx; sqrt(alphaIso)*Dy];
gtot = [g;zeros(N*N,1);zeros(N*N,1)];
tic;
fra = pcg(@(x) JTJH(x, A, -Lapl, alphaIso),A'*g,1e-6,100);
toc;

fraim = reshape(fra,N,N);
figure(5);subplot(2,2,4);imagesc(fraim); colormap(gray);title(['recon with 1-Tikononv \alpha= ',num2str(alphaIso)]);

gcorrect = zeros(185,180);
gcorrect(1:180,1:180) = fraim;
n = 128;
ig = iradon(gcorrect,angles,1,n);
figure(10);
imagesc(ig);colormap gray;title('filtered back-projection undersampled angles' );

%% repeat with edge weighted Laplacian
alphaAL = 0.01e1;
gsq = imx.^2 + imy.^2;
gmax = max(max(gsq));
T = gmax/100;
K = spdiags(exp(-reshape(gsq,[],1)/T),0,N*N,N*N);
figure(2);
subplot(2,2,2);imagesc(sqrt(gsq));colormap(gray);title('image gradient');
subplot(2,2,3);imagesc(reshape(diag(K),N,N));colormap(gray);title('\kappa');
DKD = -(Dx'*K*Dx + Dy'*K*Dy);
tic;
fra2 = pcg(@(x) JTJH(x, A, -DKD, alphaAL),A'*g,1e-6,100);
toc;

fra2im = reshape(fra2,N,N);
figure(5);subplot(2,2,3);imagesc(fraim); colormap(gray);title(['recon with 1-Tikononv \alpha= ',num2str(alphaIso)]);
figure(5);subplot(2,2,4);imagesc(fra2im); colormap(gray);title(['recon with edge-weighted Laplacian \alpha= ',num2str(alphaAL)]);

gcorrect = zeros(185,180);
gcorrect(1:180,1:180) = fra2im;
n = 128;
ig2 = iradon(gcorrect,angles,1,n);
figure(11);
subplot(2,2,1);imagesc(im);colormap gray;title('orginal image' );
subplot(2,2,2);imagesc(ig0);colormap gray;title('istotropic' );
subplot(2,2,3);imagesc(ig);colormap gray;title('diffusion process' );
subplot(2,2,4);imagesc(ig2);colormap gray;title('TV' );


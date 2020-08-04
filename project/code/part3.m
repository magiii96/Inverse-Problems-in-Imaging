%3
close all;
load('SLphan.mat');
ftrue = SLphan;
%ftrue = double(ftrue);
n = size(ftrue,1);
%3.1 Zero-order regularisation
%parameters
row = size(ftrue,1);
col = size(ftrue,2);
sigma = 0.01;
OtherParameters = struct;
OtherParameters.tol = 1e-6;
OtherParameters.maxIter = 20;
alpha = 1e2*sigma;

img2vec = @(img) reshape(img,[],1);
vec2img = @(vector) reshape(vector,row,col);

%3.1a Full range but small number of angles(zero-order)

angles2 = [0:4:179];

ATA = @(f,alpha) img2vec(iradon(radon(vec2img(f),angles2),angles2,'linear','none',1,n)) + alpha*f(:);

g = radon(ftrue,angles2);
g = g+sigma*randn(size(g));

RHS = reshape(iradon(g,angles2,'linear','none',1,n),[],1);
x0 = vec2img(pcg(@(f) ATA(f,alpha),RHS,OtherParameters.tol,OtherParameters.maxIter));
reconFBP = iradon(g,angles2,1,n); 

%Full range but small number of angles first-order Tikhonov regularisation

Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
DxT = @(x) [-x(:,1),-diff(x(:,1:end-1),1,2),x(:,end-1)];
DyT = @(x) [-x(1,:);-diff(x(1:end-1,:),1,1);x(end-1,:)];

ATA = @(f,alpha) img2vec(iradon(radon(vec2img(f),angles2),angles2,'linear', ...
    'none',1,n)) + alpha * img2vec(DxT(Dx(vec2img(f))) + ...
    DyT(Dy(vec2img(f))));

g = radon(ftrue,angles2);
g = g+sigma*randn(size(g));

RHS = reshape(iradon(g,angles2,'linear','none',1,n),[],1);
x1 = vec2img(pcg(@(f) ATA(f,alpha),RHS,OtherParameters.tol,OtherParameters.maxIter));


figure(1); clf; hold on;
subplot(2,3,1); imagesc(ftrue);title('Original Image');colorbar;
subplot(2,3,2); imagesc(g);title('Sinogram space');colorbar;
subplot(2,3,3); imagesc(reconFBP);title('Filtered Back Projection');colorbar;
subplot(2,3,4); imagesc(x0); title(['Regularised TK0, \alpha=',num2str(alpha)]);colorbar;
subplot(2,3,5); imagesc(x1); title(['Regularised TK1, \alpha=',num2str(alpha)]);colorbar;


%3.2 limited angles(zero-order)

angles3 = [0:1:44];

ATA = @(f,alpha) img2vec(iradon(radon(vec2img(f),angles3),angles3,'linear','none',1,n)) + alpha*f(:);

g = radon(ftrue,angles3);
g = g+sigma*randn(size(g));

RHS = reshape(iradon(g,angles3,'linear','none',1,n),[],1);
x0 = vec2img(pcg(@(f) ATA(f,alpha),RHS,OtherParameters.tol,OtherParameters.maxIter));
reconFBP = iradon(g,angles3,1,n); 

%Full range but small number of angles first-order Tikhonov regularisation

Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
DxT = @(x) [-x(:,1),-diff(x(:,1:end-1),1,2),x(:,end-1)];
DyT = @(x) [-x(1,:);-diff(x(1:end-1,:),1,1);x(end-1,:)];

ATA = @(f,alpha) img2vec(iradon(radon(vec2img(f),angles3),angles3,'linear', ...
    'none',1,n)) + alpha * img2vec(DxT(Dx(vec2img(f))) + ...
    DyT(Dy(vec2img(f))));

g = radon(ftrue,angles3);
g = g+sigma*randn(size(g));

RHS = reshape(iradon(g,angles3,'linear','none',1,n),[],1);
x1 = vec2img(pcg(@(f) ATA(f,alpha),RHS,OtherParameters.tol,OtherParameters.maxIter));


figure(2); clf; hold on;
subplot(2,3,1); imagesc(ftrue);title('Original Image');colorbar;
subplot(2,3,2); imagesc(g);title('Sinogram space');colorbar;
subplot(2,3,3); imagesc(reconFBP);title('Filtered Back Projection');colorbar;
subplot(2,3,4); imagesc(x0); title(['Regularised TK0, \alpha=',num2str(alpha)]);colorbar;
subplot(2,3,5); imagesc(x1); title(['Regularised TK1, \alpha=',num2str(alpha)]);colorbar;












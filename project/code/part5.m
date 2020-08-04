%5.Iterative soft-thresholding for X-ray tomography
close all;clear all;


load('SLphan.mat');
ftrue = SLphan;
n = size(ftrue,1);
row = size(ftrue,1);
col = size(ftrue,2);

% Parameters
OtherParameters = struct;
OtherParameters.tol = 10^-8;
OtherParameters.maxIter = 256;

% Function Handlers
Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
DxT = @(x) [-x(:,1),-diff(x(:,1:end-1),1,2),x(:,end-1)];
DyT = @(x) [-x(1,:);-diff(x(1:end-1,:),1,1);x(end-1,:)];

img2vec = @(img) reshape(img,[],1);
vec2img = @(vector) reshape(vector,row,col);

%limited angles 

for i=1:2
    sigma = 1;
    if i == 1
       angles = [0:4:179];
    else 
       angles = [0:1:44];
    end
    g = radon(ftrue,angles);
    g = g+sigma*randn(size(g));

    bg = iradon(g,angles,1,n); 

    alpha = 1;
    lambda = 0.0008; % step size
    ATA = @(f,alpha) img2vec(iradon(radon(vec2img(f),angles),angles,'linear','none',1,n)) + alpha*f(:);
    RHS = reshape(iradon(g,angles,'linear','none',1,n),[],1);

    % First Order Tikhonov
    f0 = pcg(@(f) ATA(f,alpha),RHS,OtherParameters.tol,OtherParameters.maxIter);

    maxIter = 256;
    tol = 0.03;

    frecon = Ite(f0,bg,ftrue,maxIter,tol,lambda,alpha,img2vec,vec2img,ATA,RHS);

    figure(4);
    subplot(3,2,i);imagesc(bg);axis off;colormap gray;...
    title(['Filter BP Noise: ', num2str(i)]);
subplot(3,2,2+i);imagesc(frecon);...
    axis off;colormap gray;title('iterate reconstruction');
subplot(3,2,2*2+i);imagesc(frecon-ftrue);...
    axis off;colormap gray;title('difference');
end

%low number of angles
for i=1:2
    sigma = 1;
    if i == 1
       angles = [0:1:179];
    else 
       angles = [0:4:179];
    end
    g = radon(ftrue,angles);
    g = g+sigma*randn(size(g));

    bg = iradon(g,angles,1,n); 

    alpha = 1;
    lambda = 0.0008; % step size
    ATA = @(f,alpha) img2vec(iradon(radon(vec2img(f),angles),angles,'linear','none',1,n)) + alpha*f(:);
    RHS = reshape(iradon(g,angles,'linear','none',1,n),[],1);

    % First Order Tikhonov
    f0 = pcg(@(f) ATA(f,alpha),RHS,OtherParameters.tol,OtherParameters.maxIter);

    maxIter = 256;
    tol = 0.03;

    frecon = Ite(f0,bg,ftrue,maxIter,tol,lambda,alpha,img2vec,vec2img,ATA,RHS);

    figure(3);
    subplot(3,2,i);imagesc(bg);axis off;colormap gray;...
    title(['Filter BP Noise: ', num2str(i)]);
subplot(3,2,2+i);imagesc(frecon);...
    axis off;colormap gray;title('iterate reconstruction');
subplot(3,2,2*2+i);imagesc(frecon-ftrue);...
    axis off;colormap gray;title('difference');

  
end



%Varying noise level
sigmas = [0.1 1 3];
angles = [0:1:179];
for i=1:size(sigmas,2)
    sigma = sigmas(i);
    g = radon(ftrue,angles);
    g = g+sigma*randn(size(g));

    bg = iradon(g,angles,1,n); 

    alpha = 1;
    lambda = 0.0008; % step size
    ATA = @(f,alpha) img2vec(iradon(radon(vec2img(f),angles),angles,'linear','none',1,n)) + alpha*f(:);
    RHS = reshape(iradon(g,angles,'linear','none',1,n),[],1);

    % First Order Tikhonov
    f0 = pcg(@(f) ATA(f,alpha),RHS,OtherParameters.tol,OtherParameters.maxIter);

    maxIter = 256;
    tol = 0.03;

    frecon = Ite(f0,bg,ftrue,maxIter,tol,lambda,alpha,img2vec,vec2img,ATA,RHS);

    figure(3);
    subplot(3,size(sigmas,2),i);imagesc(bg);axis off;colormap gray;...
    title(['Filter BP Noise: ', num2str(sigmas(i))]);
subplot(3,size(sigmas,2),size(sigmas,2)+i);imagesc(frecon);...
    axis off;colormap gray;title('iterate reconstruction');
subplot(3,size(sigmas,2),2*size(sigmas,2)+i);imagesc(frecon-ftrue);...
    axis off;colormap gray;title('difference');

end


%%CW2
%Ex1 Convolution and deconvolution 
clear all;clc;close all;
%a.
%load image
img = imread('Cameraman256.png');
%convert to float
douimg = double(img);
%Normal
imgf = douimg/max(douimg(:));
colormap(gray); imagesc(imgf); title('original image');

[row,col] = size(img);
img2vec = @(image) reshape(image,[],1);
vec2img = @(vector) reshape(vector,row,col);

%b.output the blur image Af
OtherParameters = struct;
OtherParameters.size = 5; % kernel size 
OtherParameters.sigma = 2;
BlurredIm = imblur(imgf,OtherParameters);

%add nosie 
theta = 0.1;
BlurredImnoi = BlurredIm + theta*randn(size(BlurredIm));

% plot the setup
figure();

%subplot(2,2,[1,2]); colormap(gray); imagesc(imgf); title('original image');
subplot(1,2,1); colormap(gray); imagesc(BlurredIm); title('blur image');
subplot(1,2,2); colormap(gray); imagesc(BlurredImnoi); title('blur image with noisy');


%c.

% regularization parameter alpha
alpha = 1;


% fucntion handle for (A^T A + alpha * I) * f

ATg = img2vec(imblur(BlurredImnoi,OtherParameters));

% call pcg and measure the time by tic . toc
[fa,flag,relres,iter,resvec] = pcg(@(x)ATA(x,alpha,row,col,OtherParameters),ATg);
fa = vec2img(fa);

% plot the solution
figure('Name','PCG')
subplot(2,2,1); colormap(gray); imagesc(imgf); title('original image');
subplot(2,2,2); colormap(gray); imagesc(BlurredImnoi); title('blurred image with noisy');
subplot(2,2,3); colormap(gray); imagesc(fa); title('pcg result');
subplot(2,2,4); semilogy(resvec,'g'); title('residual from Matlab PCG');
drawnow();


%d.lsqr

gv = reshape(BlurredImnoi,[],1);
g0 = [gv;zeros(row*col,1)];
AaugHandle = @(f,transposeFlag) Aaug(f,alpha,row,col,OtherParameters,transposeFlag);

[flsqr,flag,relres,iter2,resvec]= lsqr(AaugHandle,g0);
flsqr =  vec2img(flsqr);

figure('Name','LSQR')
subplot(2,2,1); colormap(gray); imagesc(imgf); title('original image');
subplot(2,2,2); colormap(gray); imagesc(BlurredImnoi); title('blurred image with noisy');
subplot(2,2,3); colormap(gray); imagesc(flsqr); title('lsqr result');
subplot(2,2,4); semilogy(resvec,'g'); title('residual from Matlab Lsqr');
drawnow();


%%EX2.choose a regularisation parameter alpha
%i discrepency principle

alphafa = @(alpha) pcgRegularization(alpha,OtherParameters,row,col,BlurredImnoi);
%ralpha2 = @(alpha) norm(imblur(alphafa(alpha),OtherParameters) -BlurredImnoi)^2;

ralpha2 = @(alpha) norm(img2vec(imblur(alphafa(alpha),OtherParameters) -BlurredImnoi))^2;
dp = @(alpha) 1/(row*col) * ralpha2(alpha) - theta^2;

logalpha = [-2:0.1:0];
alphad =  10.^logalpha;
discrepencyp = zeros(size(alphad));

for k = 1:length(alphad)
    discrepencyp(k) = dp(alphad(k));
end

figure();
semilogx(alphad,discrepencyp); title('Discrepancy Principle value');
%hold on; semilogx(alpha,log(edp),'ko');

fzeroOptions = [];
fzeroOptions.Display = 'iter';
fzeroOptions.TolX = 10^-6;
aDP = fzero(dp,10.^[-2,0],fzeroOptions);

figure('Name','alpha discrepancy');
frec = alphafa(aDP);
alphaDP = num2str(aDP);

subplot(2,3,1); colormap(gray); imagesc(imgf); title('original f_{true}');
subplot(2,3,2); colormap(gray); imagesc(BlurredImnoi); title('blur image with noisy');
subplot(2,3,3); colormap(gray); imagesc(alphafa(aDP) -imgf);title('final estimate error (f_{true}-f_{recon})'); colorbar;
subplot(2,3,4); colormap(gray); imagesc(alphafa(aDP)); title({'final reconstructed f_{recon}'; ['DP - \alpha=',num2str(aDP)]});
subplot(2,3,5); colormap(gray); imagesc(BlurredImnoi-alphafa(aDP));title('final residual (g-f_{recon})');colorbar;
subplot(2,3,6); histogram(BlurredImnoi(:)-frec(:));
drawnow();

%%ii l-curve 
millerH = @(f) 1/(2*theta^2) * norm(img2vec(imblur(f,OtherParameters)-BlurredImnoi))^2 - 0.5*alpha * norm(img2vec(f))^2;
miller    = @(alpha) millerH(alphafa(alpha));

% we first plot the Miller criterion and the L-curve for a range of alphas
alpha_vec = 10.^(-3:0.1:0);
miller_vec = zeros(size(alpha_vec));
em = miller_vec;
for k=1:length(alpha_vec)
    miller_vec(k) = miller(alpha_vec(k));
    error = img2vec(alphafa(alpha_vec(k)) -imgf);
    em(k) = norm(error)^2;

end
%Plot L-curve 

for k=1:length(alpha_vec)
    fk = alphafa(alpha_vec(k));
    llhd(k) = 0.5 * norm(img2vec(imblur(fk,OtherParameters)-BlurredImnoi))^2;
    pr(k) = 0.5*norm(img2vec(fk))^2;
end

figure('Name','Miller Criterion')
semilogx(alpha_vec,miller_vec);
hold on; semilogx(alpha_vec,em,'ko');

figure('Name','L-curve')
loglog(pr,llhd,'k'); title('L-curve');

alphaMiller = fzero(miller,10.^[-4,0],fzeroOptions);

frecL = alphafa(alphaMiller);
figure('Name','alpha miller')
subplot(2,3,1); colormap(gray); imagesc(imgf); title('original f_{true}');
subplot(2,3,2); colormap(gray); imagesc(BlurredImnoi); title('blur image with noisy');
subplot(2,3,3); colormap(gray); imagesc(alphafa(alphaMiller) -imgf);title('final estimate error (f_{true}-f_{recon})'); colorbar;
subplot(2,3,4); colormap(gray); imagesc(alphafa(alphaMiller)); title({'final reconstructed f_{recon}'; ['Miller - \alpha=',num2str(alphaMiller)]});
subplot(2,3,5); colormap(gray); imagesc(BlurredImnoi-alphafa(alphaMiller));title('final residual (g-f_{recon})');colorbar;
subplot(2,3,6); histogram(BlurredImnoi(:)-frecL(:));
drawnow();


%%Ex3 Using a regularisation term based on the spatial derivative 

gradfa = @(alpha) pcgRegularizationGradient(alpha,OtherParameters,row,col,BlurredImnoi);
ralpha2 = @(alpha) norm(img2vec(imblur(gradfa(alpha),OtherParameters) -BlurredImnoi))^2;
dp = @(alpha) 1/(row*col) * ralpha2(alpha) - theta^2;

fzeroOptions = [];
fzeroOptions.Display = 'iter';
fzeroOptions.TolX = 10^-6;
aDPGrad = fzero(dp,10.^[-1,2]);


%subplot(1,3,1); colormap(gray); imagesc(img); title('original image');
%subplot(1,3,2); colormap(gray); imagesc(BlurredImnoi); title('blur image with noisy');
%subplot(1,3,3); colormap(gray); imagesc(alphafa(aDP)); title('reconstruction, clean');

figure('Name','First Order PCG');
frecL = gradfa(aDPGrad);
subplot(2,3,1); colormap(gray); imagesc(imgf); title('original f_{true}');
subplot(2,3,2); colormap(gray); imagesc(BlurredImnoi); title('blur image with noisy');
subplot(2,3,3); colormap(gray); imagesc(gradfa(aDPGrad) -imgf);title('final estimate error (f_{true}-f_{recon})'); colorbar;
subplot(2,3,4); colormap(gray); imagesc(gradfa(aDPGrad)); title({'final reconstructed f_{recon}'; ['DP - \alpha=',num2str(aDPGrad)]});
subplot(2,3,5); colormap(gray); imagesc(BlurredImnoi-gradfa(aDPGrad));title('final residual (g-f_{recon})');colorbar;
subplot(2,3,6); histogram(BlurredImnoi(:)-frecL(:));

drawnow();

%lsqr
alphalsqr = @(alpha) lsqrRegularizationGradient(alpha,OtherParameters,row,col,BlurredImnoi);
ralphalsqr = @(alpha) norm(img2vec(imblur(alphalsqr(alpha),OtherParameters) -BlurredImnoi))^2;
lsqrdp = @(alpha) 1/(row*col) * ralphalsqr(alpha) - theta^2;

fzeroOptions = [];
fzeroOptions.Display = 'iter';
fzeroOptions.TolX = 10^-6;
alsqrDP = fzero(lsqrdp,10.^[-2,2]);

figure('Name','First Order LSQR');
%subplot(1,3,1); colormap(gray); imagesc(img); title('original image');
%subplot(1,3,2); colormap(gray); imagesc(BlurredImnoi); title('blur image with noisy');
%subplot(1,3,3); colormap(gray); imagesc(alphalsqr(alsqrDP)); title('lsqr');

frecL = alphalsqr(alsqrDP);
subplot(2,3,1); colormap(gray); imagesc(imgf); title('original f_{true}');
subplot(2,3,2); colormap(gray); imagesc(BlurredImnoi); title('blur image with noisy');
subplot(2,3,3); colormap(gray); imagesc(alphalsqr(alsqrDP) -imgf);title('final estimate error (f_{true}-f_{recon})'); colorbar;
subplot(2,3,4); colormap(gray); imagesc(alphalsqr(alsqrDP)); title({'final reconstructed f_{recon}'; ['DP - \alpha=',num2str(alsqrDP)]});
subplot(2,3,5); colormap(gray); imagesc(BlurredImnoi-alphalsqr(alsqrDP));title('final residual (g-f_{recon})');colorbar;
subplot(2,3,6); histogram(BlurredImnoi(:)-frecL(:));

drawnow();

%4.Construct an anisotropic derivative filter
% define discrete derivatives in x and y direction
Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
% compute the norm of the gradient
Normgrad= sqrt(Dx(BlurredImnoi).^2+Dy(BlurredImnoi).^2);
T = max(Normgrad(:));
gamma = exp(- Normgrad/T);

gMatrix = spdiags(gamma(:),0,row*col,row*col);
showmatrix = gMatrix(10000:12000,10000:12000);
figure('Name','filterMatrix');
imagesc(showmatrix); title('filterMatrix');

%To display 
close all
%Ex5.Iterative deblurring

alpha = 1;
fi = BlurredImnoi;
nIter = 5;
discrepancy_vec = zeros(nIter,1);

figure('Name','iterative anisotropic result')
for i = 1:nIter
    gradientNorm = sqrt(Dx(fi).^2+Dy(fi).^2);
    T = 1/alpha * max(gradientNorm(:));
    fi = AnisotropicFiltering(BlurredImnoi,OtherParameters,fi,T,alpha);
    discrepancy_vec(i) = 1/(row*col) * norm(img2vec(fi-BlurredImnoi))^2 - theta^2;
    gamma = exp(- gradientNorm/T);
    
    % plot the current iterate and its gradient image
    subplot(5,nIter,i); colormap(gray); imagesc(fi); title(['f. iter: ' int2str(i)]);
    subplot(5,nIter,i+nIter); colormap(gray); imagesc(gradientNorm); title(['grad(f), iter: ' int2str(i)]);
    subplot(5,nIter,i+2*nIter); colormap(gray); imagesc(gamma); title(['kappa. iter: ' int2str(i)]);
    subplot(5,nIter,i+3*nIter); colormap(gray); imagesc(fi-imgf); title(['residul. iter: ' int2str(i)]);
    subplot(5,nIter,i+4*nIter); histogram(BlurredImnoi(:)-fi(:));
    %subplot(6,nIter,i+5*nIter); histogram(BlurredImnoi(:)-fi(:));
    drawnow();
end



figure('Name',' Iterative deblurring');
plot(1:nIter,discrepancy_vec)

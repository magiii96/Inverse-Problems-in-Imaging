%4.A Haar wavelet denoiser 
img = im2double(imread('Cameraman256.png'));
%add noise 
sigma = 0.15;
f = img + sigma*rand(size(img));

%calculate the Haar wavelet transform and coefficients
[a,h,v,d] = haart2(f);
%plot some of coefficients 
figure(1);
subplot(3,3,1);imagesc(h{1});axis off;colormap gray;title('horizontal1');
subplot(3,3,2);imagesc(h{3});axis off;colormap gray;title('horizontal3');
subplot(3,3,3);imagesc(h{5});axis off;colormap gray;title('horizontal5');

subplot(3,3,4);imagesc(v{1});axis off;colormap gray;title('vertical1');
subplot(3,3,5);imagesc(v{3});axis off;colormap gray;title('vertical3');
subplot(3,3,6);imagesc(v{5});axis off;colormap gray;title('vertical5');

subplot(3,3,7);imagesc(d{1});axis off;colormap gray;title('diagonal1');
subplot(3,3,8);imagesc(d{3});axis off;colormap gray;title('diagonal3');
subplot(3,3,9);imagesc(d{5});axis off;colormap gray;title('diagonal5');


%reconstruct the image from the coefficients 
f_rec = ihaart2(a,h,v,d);


%compare reconstructed with original
figure(2);
subplot(1,3,1);imshow(f);title('original noisy image');
subplot(1,3,2);imshow(f_rec);title('reconstructed noisy image');
subplot(1,3,3);imshow(f_rec-f);title('difference');

%denoising for changing range.
maxRange = [1 3 5 7 8];

figure(3);
subplot(2,size(maxRange,2)+1,1);imagesc(img);...
    axis off;colormap gray;title('original');
subplot(2,size(maxRange,2)+1,1+size(maxRange,2)+1);imagesc(f);...
    axis off;colormap gray;title('blurred image');

for i=1:size(maxRange,2)
range = 1:maxRange(i);
tVal = sigma;
[hT,vT,dT] = thresholdFunction(h,v,d,range,tVal);
fdenoise = ihaart2(a,hT,vT,dT);
figure(3);
subplot(2,size(maxRange,2)+1,i+1);imagesc(fdenoise);axis off;colormap gray;...
    title(['Range: ', num2str(maxRange(i))]);
subplot(2,size(maxRange,2)+1,size(maxRange,2)+i+1+1);imagesc(f - fdenoise);...
    axis off;colormap gray;title('difference');
end

% Different tVals
tVal = sigma*[0.1 0.5 1 2 4];

figure(4);
subplot(2,size(tVal,2)+1,1);imagesc(img);...
    axis off;colormap gray;title('original');
subplot(2,size(tVal,2)+1,1+size(maxRange,2)+1);imagesc(f);...
    axis off;colormap gray;title('blurred image');


for i=1:size(tVal,2)
range = 1:size(h,2);
[hT,vT,dT] = thresholdFunction(h,v,d,range,tVal(i));
fdenoise = ihaart2(a,hT,vT,dT);
subplot(2,size(tVal,2)+1,i+1);imagesc(f_denoise);axis off;colormap gray;...
    title(['tVal: ', num2str(tVal(i))]);
subplot(2,size(tVal,2)+1,size(tVal,2) + i+1+1);imagesc(f - fdenoise);...
    axis off;colormap gray;title('difference');
end

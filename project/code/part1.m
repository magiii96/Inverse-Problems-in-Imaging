 %%Inverse project 
close all;
%1 back-projection
%1.a load image and display
load('SLphan.mat');
ftrue = SLphan;
figure(1);
imagesc(ftrue);colormap gray; title('ftrue' );
%imshow(uint8(f));
%1.b generate the random transform 
angles = [0:179];
g = radon(ftrue,angles);
figure(2);
imagesc(g);colormap gray; title('radon transform' );
%1.c unfiltered back-projection 
ori_bg = iradon(g,angles,'linear','none');
figure(3);
imagesc(ori_bg);colormap gray;title('unfiltered back-projection 1' );
%control the size of back-projection 
n = size(ftrue,1);
bg = iradon(g,angles,'linear','none',1,n);
figure(4);
imagesc(bg);colormap gray;title('unfiltered back-projection 2' );

%1.d filtered back-projection
ori_ig = iradon(g,angles);
figure(5);
imagesc(ori_ig);colormap gray; title('filtered back-projection 1' );
ig = iradon(g,angles,1,n);
figure(6);
imagesc(ig);colormap gray; title('filtered back-projection 2' );
%1.e add noise to g
sigmas = [0 0.01 0.1 0.5 1];

figure(7);
for i=1:size(sigmas,2)
    gnoise = g + sigmas(i)*randn(size(g));
    %unfiltered back projection
    bg_noise = iradon(gnoise,angles,'linear','none',1,n);
    %filtered back projection
    %bg_noise = iradon(gnoise,angles,'linear',1,n);
    subplot(1,size(sigmas,2),i);imagesc(bg_noise);...
        title(['sigma: ' num2str(sigmas(i))]);colormap gray;axis off;
end

figure(8);
for i=1:size(sigmas,2)
    gnoise = g + sigmas(i)*randn(size(g));
    ig_noise = iradon(gnoise,angles,1,n); 
    subplot(1,size(sigmas,2),i);imagesc(ig_noise);...
        title(['sigma: ' num2str(sigmas(i))]);colormap gray;axis off;
end


figure(9);
subplot(1,3,1);imagesc(A1);axis off;colormap gray;title('Angle Range: 0-180, projections: 45');













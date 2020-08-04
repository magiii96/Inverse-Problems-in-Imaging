%2.
close all;
load('SLphan.mat');
ftrue = SLphan;
n = size(ftrue,1);
%2.1 range is [0,180]. but varying number of projections
angles2 = [0:4:179];
g2 = radon(ftrue,angles2);
[numsamples,numAngles] = size(g2);

% get number of samples to create A matrix
A1 = zeros(numAngles*numsamples,n*n);

col = 1;
for i=1:n
    for j=1:n
        img = zeros(n,n);
        img(j,i) = 1;
        % take radon transform of this image
        img = radon(img,angles2);
        % reshape to column vector
        img = reshape(img,[],1);
        A1(1:size(img,1),col) = img;
        col  = col + 1;
    end
end

svd1 = svds(A1);
figure(1);
subplot(1,2,1);imagesc(A1);axis off;colormap gray;title('Angle Range: 0-180, projections: 45');


%2.2 keep the number of projectin but the range is [0,45].
angles3 = [0:1:44];
g3 = radon(ftrue,angles3);
[numsamples3,numAngles3] = size(g3);

% get number of samples to create A matrix
A2 = zeros(numAngles3*numsamples3,n*n);

col = 1;
for i=1:n
    for j=1:n
        img = zeros(n,n);
        img(j,i) = 1;
        % take radon transform of this image
        img = radon(img,angles3);
        % reshape to column vector
        img = reshape(img,[],1);
        A2(1:size(img,1),col) = img;
        col  = col + 1;
    end
end

svd2 = svds(A2);
figure(1);
subplot(1,2,2);imagesc(A2);axis off;colormap gray;title('Angle Range: 0-45, projections: 45');

figure(2);
plot(1:size(svd1,1),svd1,':*',1:size(svd2,1),svd2,':*');
legend('Angle Range: 0-180','Angle Range: 0-45');
xticks([1 2 3 4 5 6]);


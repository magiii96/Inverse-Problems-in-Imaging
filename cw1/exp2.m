%2. Singular Value Decomposition 
clear all;close all;clc
%2.a 
%set a spatial grid on the interval [−1, 1] in n equally space.
n = 20;
x2 = linspace(-1,1,n);


%2.b 
sigma = 0.2;
mu = 0;
g = 2/(sqrt(2*pi)*sigma*(n-1))*exp(-(x2 - mu).^2/(2*sigma^2));
figure1 = figure();
plot(x2,g);
saveas(figure1,'2b.png','png')

%2.c

A2 = zeros(n);
for i = 1:n
    for j = 1:n
        A2(i,j) = 2/(sqrt(2*pi)*sigma*(n-1))*exp(-(x2(i) - x2(j)).^2/(2*sigma^2)); 
    end
end


%2.d
%using imagesec
figure2 = figure();
imagesc(A2);
saveas(figure2,'2d1.png','png');

% using imwrite 

Aimg = ceil(A2/max(A2(:))*256);
figure3 = figure();
colorMap = parula(256);
imwrite(Aimg,colorMap,'2d2.png');
imshow('2d2.png');

%2.e
[U,W,V] = svd(A2);

%check A B are equal matrix.
B2 = U*W*V';
check = norm(A2-B2);

%2.f
%using sparse to construct the pseudoinverse of W.
winv = spdiags(1./diag(W),0,n,n);
figure4 = figure(); 
spy(winv);
saveas(figure4,'2f1_20.png','png');

%check W*invw = n * n Identity matrix
M1 = W*winv;

%check the pseudoinverse of A.
C2 = V*winv*U';
Ainv = pinv(A2);
M2 = C2*A2;
%should be zero when n=10,20, but not zeron when n =100.
check2 = norm(Ainv - C2);

%plot diag(W) on a logarithmic scale on the y-axis
figure5 = figure();
semilogy(diag(W));
saveas(figure5,'2f2_20.png','png');

%2.g
%repeat last two steps when n = 20,100;
%plot the first 9 columns of V , the last 9 columns of V
figure6 = figure();
subplot(1,2,1)
plot(V(:,1:9))
subplot(1,2,2)
plot(V(:,end-9:end));
saveas(figure6,'2g.png','png');


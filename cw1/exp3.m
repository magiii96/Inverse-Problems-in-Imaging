clear all;close all;clc
%3.a
n3 = 100;
x3 = linspace(-1,1,n3);
xa = zeros(1,n3);

for a = 1:n3
    if (-0.95 < x3(a)) && (x3(a)<= -0.6)
       xa(a) = 1;
    end
end

xb = zeros(1,n3);
for b = 1:n3
    if (-0.6 < x3(b)) && (x3(b)<= -0.2)
       xb(b) = 1;
    end
end

xc = zeros(1,n3);
for c = 1:n3
    if (-0.2 < x3(c)) && (x3(c) <= 0.2)
       xc(c) = 1;
    end
end

xd = zeros(1,n3);
for d = 1:n3
    if (0.4 < x3(d)) && (x3(d) <= 0.6)
       xd(d) = 1;
    end
end

xe = zeros(1,n3);
for e = 1:n3
    if (0.6 < x3(e)) && (x3(e) <= 1)
       xe(e) = 1;
    end
end

f = xa+0.2*xb-0.5*xc+0.7*xd-0.7*xe;
figure1 = figure();
plot(x3,f);
%saveas(figure1,'3a.png','png')

%3.b
sigma3 = 0.2;
A3 = zeros(n3,n3);
for i = 1:n3
    for j = 1:n3
        A3(i,j) = 2/(sqrt(2*pi)*sigma3*(n3-1))*exp(-(x3(i) - x3(j)).^2/(2*sigma3^2)); 
    end
end

%using imagesec
figure2 = figure();
imagesc(A3);
%saveas(figure2,'3b1_01.png','png');

[U3,W3,V3] = svd(A3);


%plot diag(W) on a logarithmic scale on the y-axis
figure3 = figure();
semilogy(diag(W3));
%saveas(figure3,'3b2_01.png','png');

%3.c verify following the normal
n_diag = sqrt(-log(diag(W3)));
figure4 = figure();
plot(n_diag)
saveas(figure4,'3c_02.png','png');

%3.d
ft = f';
result = A3 * ft;
resultt = result';

figure5 = figure();
plot(x3,resultt);
saveas(figure5,'3d_02.png','png');




%3.e
g = 2/(sqrt(2*pi)*sigma3*(n3-1))*exp(-(x3).^2/(2*sigma3^2));

Fx = fftshift(fft(fftshift(f)));
gx = fftshift(fft(fftshift(g)));
resultx =Fx .* gx;
resultx = fftshift(ifft(fftshift(resultx)));

figure6 = figure();
plot(x3,resultx);
saveas(figure6,'3e_02.png','png');

figure7 = figure();
plot(x3,resultx);
hold on 
%compare two convolution 
plot(x3,resultt);
legend('Fourier space','convolution matrix')
saveas(figure7,'3e_c02.png','png');



%3.f
periodicA1 = circshift(A3,[50 50]);
periodicA = zeros(n3,n3);
for i = 1:n3
    for j = 1:n3
        periodicA(i,j) = max(periodicA1(i,j),A3(i,j));
    end
end
figure8 = figure();
imagesc(periodicA);
saveas(figure8,'3f1','png');


resultp = periodicA * ft;
resultp = resultp';

figure9 = figure();
plot(x3,resultx);
hold on 
%compare two convolution 
plot(x3,resultp);
legend('Fourier space','convolution matrix with periodic boundary conditions ')
saveas(figure9,'3f_c02.png','png');











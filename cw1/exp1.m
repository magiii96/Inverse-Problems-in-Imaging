%1.a define the function phi
%x is a vector with length 2 and p is a scalar.
p = 1;
phi = @(x,p) sum(abs(x).^p);

%1.b use fmincon to solve our problem 
A = [1,2];
b = 5;


x = fmincon(@(x)phi(x,p),[0,0],[],[],A,b);

%1.c plot our solution
px = linspace(-0.5,2);
py = (-1/2)*(px-5);
%lpy = linspace(1,3);
pointx = [0 5/9 1 1.197661529219968 1.306019319782959 1.373998402776017 1.420517627418583 ];
pointy = [2.5 20/9 2 1.901168929987217 1.846990337069687 1.813000795573158 1.789741153594409];
plot(px,py,pointx,pointy,'*');


%1.d compute the Moore-Penros inverse.
Ai = pinv(A);
Ainv = A.'*(inv( A * A.'));
%we can see that we get the same inverse matrix
xmp = Ainv*b;

%we will get different solution 
xb = A\b;
plot(px,py,pointx,pointy,'*',xmp(1),xmp(2),'b--o');

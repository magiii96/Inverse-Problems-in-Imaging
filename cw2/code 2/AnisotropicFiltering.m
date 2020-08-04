function falpha = AnisotropicFiltering(g,OtherParameters,fEdges,T,alpha)
[row, col] = size(g);
%function handles
A  = @(x)  imfilter(x,Aker,'circular');
Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
DxT = @(x) [-x(:,1),-diff(x(:,1:end-1),1,2),x(:,end-1)];
DyT = @(x) [-x(1,:);-diff(x(1:end-1,:),1,1);x(end-1,:)];
im2vec = @(image) reshape(image,[],1);
vec2im = @(vector) reshape(vector,row,col);
kappa = exp(- sqrt(Dx(fEdges).^2+Dy(fEdges).^2)/T);

% reshape f to an image format
ATAalphaGradKappa = @(f) im2vec(imblur(imblur(vec2im(f),OtherParameters),OtherParameters)...
+alpha * (DxT(kappa .* Dx(vec2im(f))) + DyT(kappa.* Dy(vec2im(f)))));
rightHandSide = reshape(imblur(g,OtherParameters),[],1);

[falpha,flag] = pcg(ATAalphaGradKappa,rightHandSide);
falpha = reshape(falpha,row,col);
end
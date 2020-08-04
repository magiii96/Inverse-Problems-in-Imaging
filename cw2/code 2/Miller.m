function [mil,e] = Miller(alpha,OtherParameters,row,col,BlurredImnoi,g,A,ftrue,sigma,RH)

% calculate discrepency value
ATAresult = @(x) ATA(x,alpha,row,col,OtherParameters);
ATg = reshape(imblur(BlurredImnoi,OtherParameters),[],1);
%tol = [];
%maxIter = 64;
fa = pcg(ATAresult,ATg);
fa = reshape(fa,[row,col]);


mil = 1/(2*theta^2) * norm(img2vec(imblur(f,OtherParameters)-BlurredImnoi))^2 - 0.5*alpha * norm(img2vec(f))^2;

n = size(A,1);
fi = (A'*A + alpha*RH)\(A'*g); % pseudo inverse solution
r = BlurredImnoi -A*fi;
dp = r'*r./n / sigma^2;
%calculate Prior;
pp = 0.5*fi'*RH*fi/n;
mil = 0.5*dp - pp;
df = fi - ftrue;
e = df'*df;

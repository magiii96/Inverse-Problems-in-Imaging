%pcg
function flsqr = lsqrRegularizationGradient(alpha,OtherParameters,row,col,BlurredImnoi)

gv = reshape(BlurredImnoi,[],1);
g0 = [gv;zeros(row*col,1)];

AaugHandle = @(f,transposeFlag) AaugGrad(f,transposeFlag,alpha,row,col,OtherParameters);


[flsqr,flag] = lsqr(AaugHandle,g0);

flsqr = reshape(flsqr,[row,col]);

end

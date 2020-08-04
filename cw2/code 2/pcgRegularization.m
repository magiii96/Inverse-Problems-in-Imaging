%pcg
function fa = pcgRegularization(alpha,OtherParameters,row,col,BlurredImnoi)

ATAresult = @(x) ATA(x,alpha,row,col,OtherParameters);

ATg = reshape(imblur(BlurredImnoi,OtherParameters),[],1);
%tol = [];
%maxIter = 100;
fa = pcg(ATAresult,ATg);

fa = reshape(fa,[row,col]);

end

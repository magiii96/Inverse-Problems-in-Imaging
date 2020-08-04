%pcg
function fa = pcgRegularizationGradient(alpha,OtherParameters,row,col,BlurredImnoi)

ATAresult = @(x) ATAgrad(x,alpha,row,col,OtherParameters);

ATg = reshape(imblur(BlurredImnoi,OtherParameters),[],1);

fa = pcg(ATAresult,ATg);

fa = reshape(fa,[row,col]);

end

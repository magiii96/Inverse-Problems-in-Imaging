
function z = ATA(f,alpha,row,col,OtherParameters)
    f = reshape(f,row,col);
    y = imblur(f,OtherParameters);
    z = imblur(y,OtherParameters) + alpha*f;
    z = z(:);
end
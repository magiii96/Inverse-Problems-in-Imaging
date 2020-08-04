function z=Aaug(f,alpha,row,col,OtherParameters,transposeFlag)
    switch transposeFlag
        case 'notransp'
        % implementation of the augmented matrix multiplication
        f = reshape(f,row,col);
        Af = imblur(f,OtherParameters);
        z = [Af(:);sqrt(alpha)*f(:)];
        case 'transp'
        % implementation of the transposed augmented matrix multiplication
        f1 = reshape(f(1:row*col),row,col);
        f2 = reshape(f((row*col+1):end),row,col);
        y = imblur(f1,OtherParameters);
        z = y(:) + sqrt(alpha)*f2(:);
        otherwise
        error('input transposeFlag has to be "transp" or "notransp"')
    end
end

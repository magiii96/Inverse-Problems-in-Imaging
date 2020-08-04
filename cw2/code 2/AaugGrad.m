function z=AaugGrad(f,transposeFlag,alpha,row,col,OtherParameters)
    Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
    Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
    DxT = @(x) [-x(:,1),-diff(x(:,1:end-1),1,2),x(:,end-1)];
    DyT = @(x) [-x(1,:);-diff(x(1:end-1,:),1,1);x(end-1,:)];
    switch transposeFlag
        case 'notransp'
        % implementation of the augmented matrix multiplication
        f = reshape(f,row,col);
        y = imblur(f,OtherParameters);
        gradf = DxT(Dx(f))+DyT(Dy(f));
        z = [y(:);sqrt(alpha)*gradf(:)];
        case 'transp'
        % implementation of the transposed augmented matrix multiplication
        f1 = reshape(f(1:row*col),row,col);
        f2 = reshape(f((row*col+1):end),row,col);
        y = imblur(f1,OtherParameters);
        gradf2 = DxT(Dx(f2))+DyT(Dy(f2));
        z = y(:) + sqrt(alpha)*gradf2(:);
        otherwise
        error('input transposeFlag has to be "transp" or "notransp"')
    end
end


function z = ATAgrad(f,alpha,row,col,OtherParameters)
    Dx = @(x)  [diff(x,1,2),zeros(size(x,1),1)];
    Dy = @(x)  [diff(x,1,1);zeros(1,size(x,2))];
    DxT = @(x) [-x(:,1),-diff(x(:,1:end-1),1,2),x(:,end-1)];
    DyT = @(x) [-x(1,:);-diff(x(1:end-1,:),1,1);x(end-1,:)];

    f = reshape(f,row,col);
    y = imblur(f,OtherParameters);
    z = imblur(y,OtherParameters) + alpha*(DxT(Dx(f))+DyT(Dy(f)));
    z = z(:);
end
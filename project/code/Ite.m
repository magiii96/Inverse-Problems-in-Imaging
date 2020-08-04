%% Functions
function fk1 = MainIteration(f0,f_bg,f,maxIter,tol,lambda,alpha,img2vec,vec2img,ATA,RHS)
tVal = alpha*lambda;
iter = 0;
% main iteration
fk = f0;
fig = figure;
while iter < maxIter
    f_temp = fk - lambda * (ATA(fk,alpha) - RHS);
    [a,h,v,d] = haart2(vec2img(f_temp));
    range = 1:size(h,2);
    [hT,vT,dT] = thresholdFunction(h,v,d,range,tVal);
    fk1 = img2vec(ihaart2(a,hT,vT,dT));
    
    % stopping condition
    if abs(norm(fk1 - fk)) <= tol
        close(fig);
        fk1 = vec2img(fk1);
        break;
    end
    disp(['Iteration: ' num2str(iter)]);
    disp(['Norm change in image: ' num2str(abs(norm(fk1 - fk)))]);
    figure(1);
    subplot(1,3,1);imagesc(f_bg);axis off;colormap gray;
    subplot(1,3,2);imagesc(vec2img(fk1));axis off;colormap gray;
    subplot(1,3,3);imagesc(vec2img(f) - vec2img(fk1));axis off;colormap gray;
    pause(0.01);
    fk = fk1;
    iter = iter + 1;
end

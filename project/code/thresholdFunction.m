%threshold function
function [hT,vT,dT] = thresholdFunction(h,v,d,range,tVal)
hT = cell(1,size(range,2));
vT = cell(1,size(range,2));
dT = cell(1,size(range,2));
for i=1:size(range,2)
    hT{i} = wthresh(h{range(i)},'s',tVal);
    vT{i} = wthresh(v{range(i)},'s',tVal);
    dT{i} = wthresh(d{range(i)},'s',tVal);
end
end
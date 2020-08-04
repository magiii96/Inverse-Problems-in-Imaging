%img blur function
function Af = imblur(f,OtherParameters)
    Aker = fspecial('gaussian',OtherParameters.size,OtherParameters.sigma);
    Af = imfilter(f, Aker,'circular');
end
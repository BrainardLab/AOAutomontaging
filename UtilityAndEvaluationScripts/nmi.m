function v = nmi(x, y)
% Compute nomalized mutual information I(x,y)/sqrt(H(x)*H(y)).
% Written by Michael Chen (sth4nth@gmail.com).
    assert(numel(x) == numel(y));
    n = numel(x);
    x = reshape(x,1,n);
    y = reshape(y,1,n);
    
    l = min(min(x),min(y));
    x = x-l+1;
    y = y-l+1;
    k = max(max(x),max(y));

    idx = 1:n;
    Mx = sparse(idx,x,1,n,k,n);
    My = sparse(idx,y,1,n,k,n);
    Pxy = nonzeros(Mx'*My/n); %joint distribution of x and y
    Hxy = -dot(Pxy,log2(Pxy+eps));

    Px = mean(Mx,1);
    Py = mean(My,1);

    % entropy of Py and Px
    Hx = -dot(Px,log2(Px+eps));
    Hy = -dot(Py,log2(Py+eps));

    % mutual information
    %check for corner cases
    if(Hx <=0 && Hy <=0 && Hxy <=0)%if both images are constant
       v = 1;
    elseif(Hx <=0 || Hy <=0 || Hxy <=0)%if one image is constant and the other is not
       v = 0;
    else
        MI = Hx + Hy - Hxy;
        % normalized mutual information
        %v = 2*MI/(Hx + Hy);
        v = sqrt((MI/Hx)*(MI/Hy)) ;
    end
end   
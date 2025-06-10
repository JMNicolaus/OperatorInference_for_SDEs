function [y, dxdt, ind] = central_finite_differences(x,h,acc,dim)
% coeficients for central FD taken from:
% https://en.m.wikipedia.org/wiki/Finite_difference_coefficient
    if ~exist("acc","var")
        acc = 2;
    end
    if ~exist("dim","var")
      dim = 2;
    end

    switch acc
        case 2
            indx = -1:1;
            coef = 1/h*[-1/2 0 1/2];
        case 4
            indx = -2:2;
            coef = 1/h*[1/12 -2/3 0 2/3 -1/12];
        case 6
            indx = -3:3;
            coef = 1/h*[-1/60 3/20 -3/4 0 3/4 -3/20	1/60 ];
        case 8
            indx = [-4:4];
            coef = 1/h*[1/280 -4/105 1/5 -4/5 0 4/5 -1/5 4/105 -1/280 ];
        otherwise
            error("acc can only be 2 4 6 or 8!")
    end
    
    switch dim 
      case 1
        coef = reshape(coef,[],1);
      case 2 
        %pass
      case 3
        coef = reshape(coef,1,1,[]);
      otherwise 
        error("dim can only be 1,2,3")
    end
    
    p = numel(coef);
    ind = 1:size(x,dim);
    ind = ind((p-1)/2+1:end-(p-1)/2);
    dxdt =[];
    y = [];

    for ii = ind
      switch dim
        case 1
         dX = sum(coef.*x(ii+indx,:,:),dim);
         Xii = x(ii,:,:);
        case 2
          dX = sum(coef.*x(:,ii+indx,:),dim);
          Xii = x(:,ii,:);
        case 3
          dX = sum(coef.*x(:,:,ii+indx),dim);
          Xii = x(:,:,ii);
      end
        dxdt = cat(dim,dxdt,dX);
        y = cat(dim,y,Xii);
    end
end
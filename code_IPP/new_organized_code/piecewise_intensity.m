% piece intensity function estimate.

function [lam_star,lambda,xarg,yarg,zarg]=piecewise_intensity(n_bins,T,S,plt,xgrid,ygrid,zgrid)
% Input:
% n_bins: the total number of bins used to do piecewise intensity
% estimation, if n_bin=[], we shall use the xgrid, ygrid and zgrid
% specified outside the function.
% T: the range of the domains, each row indicate the smallest and maximum
% bound of the dimension.
% S: the location of the events, n by dim matrix.
% plt: if 1, plot the estimated intensity.
% xgrid, ygrid, zgrid: if n_bins=[], this grid will determine the bins. 

T_dim=size(T,1);
T_area=prod(range(T,2));
if isempty(n_bins)==0
    unitbin_area=T_area/n_bins;
end
if T_dim==2    
    
   if isempty(n_bins)==0
        Tx=range(T(1,:));
        Ty=range(T(2,:));    
        ratio_yx=Ty/Tx;
        stepsize_x=sqrt(unitbin_area/ratio_yx);
        stepsize_y=ratio_yx*stepsize_x;
        xgrid=linspace(T(1,1),T(1,2),Tx/stepsize_x);
        ygrid=linspace(T(2,1),T(2,2),Ty/stepsize_y);    
   end
   
   lambda=NaN(length(xgrid)-1,length(ygrid)-1);
    for i=1:(length(xgrid)-1)
        for j=1:(length(ygrid)-1)
            id=find(((S(:,1)>=xgrid(i))&(S(:,1)<xgrid(i+1)))&((S(:,2)>=ygrid(j))&(S(:,2)<ygrid(j+1))));
            if isempty(n_bins)
               unitbin_area=(xgrid(i+1)-xgrid(i))*(ygrid(j+1)-ygrid(j));
            end
            if isempty(id)
               lambda(i,j)=0;
            else
               lambda(i,j)=length(id)/unitbin_area; 
            end        
        end
    end
    
    xarg=(xgrid(1:end-1)+xgrid(2:end))/2;
    yarg=(ygrid(1:end-1)+ygrid(2:end))/2;
    zarg=[];      
    if plt==1
        n_bins=numel(lambda);
        figure() 
          imagesc(xarg,yarg, lambda');             
          title(['nbin= ', int2str(n_bins), '   \lambda^*=',num2str(max(max(lambda)))]);
    end
    lam_star=max(max(lambda));
   
elseif T_dim==3
    
    if isempty(n_bins)==0
        Tx=range(T(1,:));
        Ty=range(T(2,:));    
        Tz=range(T(3,:));
        ratio_yx=Ty/Tx;
        ratio_zx=Tz/Tx;
        stepsize_x = nthroot(unitbin_area/ratio_yx/ratio_zx, 3);
        stepsize_y = ratio_yx*stepsize_x;
        stepsize_z = ratio_zx*stepsize_x;
        xgrid=linspace(T(1,1),T(1,2),Tx/stepsize_x);
        ygrid=linspace(T(2,1),T(2,2),Ty/stepsize_y);
        zgrid=linspace(T(3,1),T(3,2),Ty/stepsize_z);
    end
    
    lambda=NaN(length(xgrid)-1,length(ygrid)-1,length(zgrid)-1);
    for i=1:(length(xgrid)-1)
        for j=1:(length(ygrid)-1)
            for k=1:(length(zgrid)-1)
                id=find(((S(:,1)>=xgrid(i))&(S(:,1)<xgrid(i+1)))&((S(:,2)>=ygrid(j))&(S(:,2)<ygrid(j+1)))&((S(:,3)>=zgrid(k))&(S(:,3)<zgrid(k+1))));
                 if isempty(n_bins)
                    unitbin_area=(xgrid(i+1)-xgrid(i))*(ygrid(j+1)-ygrid(j))*(zgrid(k+1)-zgrid(k));
                 end
                if isempty(id)
                   lambda(i,j,k)=0;
                else
                   lambda(i,j,k)=length(id)/unitbin_area; 
                end                              
            end
       end    
    end
    
    xarg=(xgrid(1:end-1)+xgrid(2:end))/2;
    yarg=(ygrid(1:end-1)+ygrid(2:end))/2;
    zarg=(zgrid(1:end-1)+zgrid(2:end))/2;
    
    lam_star=max(max(max(lambda)));
    if plt==1
        id=find(lambda~=0);        
        dx=length(xgrid)-1;
        dy=length(ygrid)-1;
        dz=length(zgrid)-1;
                       
        figure()                                 
        for i=1:length(id)
            [idx,idy,idz]=id_to_coordinate3d(id(i),dx,dy,dz);    
            scatter3(xarg(idx),yarg(idy),zarg(idz),lambda(idx,idy,idz)*10,ceil(lambda(idx,idy,idz)/lam_star*256),'filled');            
            hold on            
        end 
        daspect([range(T(1,:)),range(T(2,:)),range(T(3,:))]);        
        hold off
    end
    
   
end




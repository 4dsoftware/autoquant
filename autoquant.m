function [pk,pkreg,flag1] = autoquant(y0,window,beta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AUTOQUANT automatically finds the peaks. Note that the unit is index.
%
% Input:
% y0 = original signal
% window = window size of SG filter
% beta = the upper threshold of noise normalized by the maximum of original
% signal, normally 0 ~ 0.1
%
% Output: 
% pk = the index of the point corresponding to peak location
% pkreg = the indexes of peak end points
% flag1 = incomplete peak flag
%
% Copyright: Di Du, Ph.D., ExxonMobil Research and Engineering
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculation of first and second derivative, was input in the old version,
%updated on 11/19/2018
order = 2; %DW statistic, normally fixed at 2
y1 = sgolayfilt(deriv(y0),order,window); %first derivative
y2 = sgolayfilt(deriv(y1),order,window); %second derivative

%initialization
flag1 = 0; %flag for incomplete peak
lengy = length(y0); %length of the input time series
thres = max(abs(y1))*beta; %upper threshold of noise in y1
valley_thre = 0.3; %used in y2
mag_limit = 10; %used in y0
peak_dist = 30; %used to exclude peaks too close to each other
cage_cap = 12; %cage capacity
rev_tol = 5; %reverse trending tolerance

%find valleys of y2
[lows,ind] = findpeaks(-y2);
pk = ind(lows>(valley_thre*max(-y2))); %count valleys deep enough

%check if the peak is incomplete
[~,miny2] = min(y2);
if isempty(find(pk == miny2, 1)) == 1
    flag1 = 1; %incomplete peak
    pk = miny2; %still continue to integration
end

%filter outstanding peaks in y0
bad_pk = [];
for i = 1:length(pk)
    if (y0(pk(i)) - min(y0)) < mag_limit
        bad_pk = [bad_pk i];        
    end
end
pk(bad_pk) = [];
if isempty(pk) == 1 
    %if all peaks are below limit, then take the most outstanding one
    pk = miny2;  
end

%peak lumping
pk = sort(pk);
lpk = length(pk);
if lpk > 1
    i = 1;
    while i < lpk 
        if (pk(i+1) - pk(i)) < peak_dist
            pk(i + heaviside(y0(pk(i))-y0(pk(i+1)))) = [];
        else
            i = i + 1;
        end
        lpk = length(pk);
    end
end
pkreg = zeros(length(pk),2);

%find peak region
for i = 1:length(pk)
    if pk(i) > 1 && pk(i) < lengy
        
    pk_left = pk(i) - heaviside(-y1(pk(i)));
    pk_right = pk(i) + heaviside(y1(pk(i)));
    
    %search left
    in_cage = 0;
    in_rev = 0;
    init = 0;
    if pk_left > 2
    for j = 1 : (pk_left - 2)
        if y1(pk_left-j) > thres && y1(pk_left-j-1) < thres
            init = 1;
        end
        if init == 1
            if y1(pk_left-j) < thres && y1(pk_left-j) > -thres %fall into cage
                in_cage = in_cage + 1;
            end
            if y1(pk_left-j) < -thres %through -thres
                in_rev = in_rev + 1;
            end
            if in_cage > cage_cap %cage is full
                pkreg(i,1) = pk_left - j;
                break;
            end                
            if in_rev > rev_tol %reverse pool is full 
                pkreg(i,1) = pk_left - j + rev_tol;
                break;
            end
        end
    end  
    if j == (pk_left - 2)
        pkreg(i,1) = 3;
    end
    else
        pkreg(i,1) = 3;
    end
    
    %search right
    in_rev = 0;
    init = 0;
    if pk_right < lengy - 1
    for j = 1 : (lengy-pk_right-1)
        if y1(pk_right+j) < -thres && y1(pk_right+j+1) > -thres
            init = 1;
        end
        if init == 1
        if y1(pk_right+j) > thres %through thres
            in_rev = in_rev + 1;
        end
        if in_rev > rev_tol %through (-thres) 
            pkreg(i,2) = pk_right + j + 1 - rev_tol;
            break;
        end
        end 
    end
    if j == (lengy-pk_right-1)
        pkreg(i,2) = lengy - 2;
    end
    else
        pkreg(i,2) = lengy - 2;
    end
    
    elseif pk(i) == 1
    pkreg(i,1) = 3;
    pk_right = pk(i) + heaviside(y1(pk(i)));
    in_rev = 0;
    init = 0;
    if pk_right < lengy - 1
    for j = 1 : (lengy-pk_right-1)
        if y1(pk_right+j) < -thres && y1(pk_right+j+1) > -thres
            init = 1;
        end
        if init == 1
        if y1(pk_right+j) > thres %through thres
            in_rev = in_rev + 1;
        end
        if in_rev > rev_tol %through (-thres) 
            pkreg(i,2) = pk_right + j + 1 - rev_tol;
            break;
        end
        end 
    end
    if j == (lengy-pk_right-1)
        pkreg(i,2) = lengy - 2;
    end
    else
        pkreg(i,2) = lengy - 2;
    end
        
    elseif pk(i) == lengy
    pkreg(i,2) = lengy - 2;
    pk_left = pk(i) - heaviside(-y1(pk(i)));
    in_cage = 0;
    in_rev = 0;
    init = 0;
    if pk_left > 2
    for j = 1 : (pk_left - 2)
        if y1(pk_left-j) > thres && y1(pk_left-j-1) < thres
            init = 1;
        end
        if init == 1
            if y1(pk_left-j) < thres && y1(pk_left-j) > -thres %fall into cage
                in_cage = in_cage + 1;
            end
            if y1(pk_left-j) < -thres %through -thres
                in_rev = in_rev + 1;
            end
            if in_cage > cage_cap %cage is full
                pkreg(i,1) = pk_left - j;
                break;
            end                
            if in_rev > rev_tol %reverse pool is full 
                pkreg(i,1) = pk_left - j + rev_tol;
                break;
            end
        end
    end  
    if j == (pk_left - 2)
        pkreg(i,1) = 3;
    end
    else
        pkreg(i,1) = 3;
    end
    
    end
end


    

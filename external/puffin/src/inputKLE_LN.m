function [ KL ] = inputKLE_LN( ndim, Kmg, Ksg )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


if ndim == 1
    mul = exp(Kmg+ (0.5*Ksg^2));  %% Mean for LogNormal RV
    str = strcat('cijk_L3_10');   %%: Loads Cijk usable for Cijk_nkle_nord
    load(str)
    
    KL(1) = mul*1;
    KL(2) = mul*Ksg;
    KL(3) = mul*(Ksg^2/2);
    
elseif ndim == 2
    
    mul = exp( Kmg+ (0.5*(Ksg^2+Ksg^2)) ); %% Mean for LogNormal RV
    
    KL(1) = mul*1;
    KL(2) = mul*Ksg;
    KL(3) = mul*Ksg;
    KL(4) = mul*(Ksg^2/2);
    KL(5) = mul*(Ksg^2/1);
    KL(6) = mul*(Ksg^2/2);
    
elseif ndim == 3
    
    mul = exp( Kmg+ (0.5*(Ksg^2+Ksg^2+Ksg^2)) ); %% Mean for LogNormal RV
    
    KL(1) = mul*1;
    KL(2) = mul*Ksg;
    KL(3) = mul*Ksg;
    KL(4) = mul*Ksg;
    KL(5) = mul*(Ksg^2/2);
    KL(6) = mul*(Ksg^2/1);
    KL(7) = mul*(Ksg^2/1);
    KL(8) = mul*(Ksg^2/2);
    KL(9) = mul*(Ksg^2/1);
    KL(10) = mul*(Ksg^2/2);
    
elseif ndim == 5
    
    mul = exp( Kmg+ (0.5*(Ksg^2+Ksg^2+Ksg^2+Ksg^2+Ksg^2)) ); %% Mean for LogNormal RV
    
    KL(1) = mul*1;
    KL(2) = mul*Ksg;
    KL(3) = mul*Ksg;
    KL(4) = mul*Ksg;
    KL(5) = mul*Ksg;
    KL(6) = mul*Ksg;
    KL(7) = mul*(Ksg^2/2);
    KL(8) = mul*(Ksg^2/1);
    KL(9) = mul*(Ksg^2/1);
    KL(10) = mul*(Ksg^2/1);
    KL(11) = mul*(Ksg^2/1);
    KL(12) = mul*(Ksg^2/2);
    KL(13) = mul*(Ksg^2/1);
    KL(14) = mul*(Ksg^2/1);
    KL(15) = mul*(Ksg^2/1);
    KL(16) = mul*(Ksg^2/2);
    KL(17) = mul*(Ksg^2/1);
    KL(18) = mul*(Ksg^2/1);
    KL(19) = mul*(Ksg^2/2);
    KL(20) = mul*(Ksg^2/1);
    KL(21) = mul*(Ksg^2/2);

end

end


function [cnfact] = Norcoeffi(M)
% cnfact  normalization coefficient for qammod.m function/ because we need
% a fixed transmitted power for different constellation size
% M is the constellation size
% cnfact is the coefficent that should be multiplied after qammod or
% divided before qamdemod
% Date:26-Oct-2020
% Cuiwei HE
if M==2
    x=[-1,1];
    y=0; 
elseif M==8
    x=-(sqrt(16)-1):2:sqrt(16)-1;
    y=[-1,1];    
else
    x=-(sqrt(M)-1):2:sqrt(M)-1;
    y=-(sqrt(M)-1):2:sqrt(M)-1;
end

sum=0;

for i=1:1:length(x)
    for j=1:1:length(y)
        sum=sum+(x(i)^2+y(j)^2);
    end
end

ave=sum/M;

cnfact=sqrt(2/ave);

end


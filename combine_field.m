clear all


imax=24;
jmax=24;
nI=3;
nJ=3;
fullfield=zeros(imax*nI,jmax*nJ);

for I=1:nI
    for J=1:nJ
        filename=['field_I',num2str(I),'_J',num2str(J),'.txt']
        subfield=load(filename);
        i=(I-1)*imax;
        j=(J-1)*jmax;
        size(fullfield(i+1:i+imax,j+1:j+jmax));
        size(subfield(:,:));
        
        fullfield(i+1:i+imax,j+1:j+jmax)=subfield(:,:);
        
    end
end

figure(3)
surf(fullfield)
size(fullfield)

save fullfield

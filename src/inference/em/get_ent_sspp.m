function ent = get_ent_sspp( stats )


covsmth = stats.covsmth;
[stadim,d2,totsamp] = size(covsmth);
sumdet = 0;
for nsamp = 1:totsamp
    sumdet = sumdet + det(covsmth(:,:,nsamp));
end
ent = 0.5*( (stadim*totsamp)*log(2*pi) + log(sumdet) + stadim*totsamp);
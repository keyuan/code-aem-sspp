function ent = get_ent_sspp_s( stats )

covsmth = stats.covsmth;
ent     = 0.5*sum(log(2*pi)+log(covsmth(2:end))+1);
Vdistort = sum(sum(map.*double(Vmask)/100))/sum(sum(map))
Ddistort = sum(sum(map.*double(Dmask)/100))/sum(sum(map))
Sdistort = sum(sum(map.*double(Smask)/100))/sum(sum(map))

fprintf('USA:\n');
Vdistort = sum(sum(map(Nmask(:)==237).*double(Vmask(Nmask(:)==237))/100))/sum(sum(map(Nmask(:)==237)))
Ddistort = sum(sum(map(Nmask(:)==237).*double(Dmask(Nmask(:)==237))/100))/sum(sum(map(Nmask(:)==237)))
Sdistort = sum(sum(map(Nmask(:)==237).*double(Smask(Nmask(:)==237))/100))/sum(sum(map(Nmask(:)==237)))
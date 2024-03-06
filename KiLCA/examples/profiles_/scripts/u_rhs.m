function dy = u_rhs(r, y)

global a alpha theta_0 mu_0;

%y = [rBth Bz p]

mu = 2.0/a * theta_0 * (1.0 - (r/a).^alpha);

for k=1:length(r),

  if r(k) > 0,

    Bth = y(1,k)./r(k);
    Bz = y(2,k);
    p = 0.0;

    B2 = Bth.^2 + Bz.^2;

    dy(3,k) = 0.0; %zero gas pressure

    dy(2,k) = - mu(k).*Bth - mu_0 * Bz./B2 * dy(3,k);    %d(rBth)

    dy(1,k) = r(k).*(mu(k).*Bz - mu_0 * Bth./B2 * dy(3,k));    %d(rBth)

  else

    dy(1:3,k) = 0.0;

  end

end

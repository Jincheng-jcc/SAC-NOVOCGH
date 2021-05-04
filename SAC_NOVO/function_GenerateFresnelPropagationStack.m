function [HStack]=function_GenerateFresnelPropagationStack(Nx,Ny,z, lambda, psx, psy, useGPU,NA,RI)
% Lambda, ps, z has unit meter.
% z is the distance from focal plane.
% lambda is wavelength
% ps is pixel size.
% RI is the refractive index of the objective
cx=(1:Nx) - (floor(Nx/2)+1);
cy=(1:Ny) - (floor(Ny/2)+1);

if useGPU
    cx = gpuArray(cx); cy = gpuArray(cy);
end

% 1.traditional transfer function in NOVO
% 1)without any modification
%  [us, vs]=meshgrid(cx, cy);
% us=us/Nx/psx; vs=vs/Ny/psy;
% HStack = exp(1i*lambda*pi*(us.^2+vs.^2)* z/2);
% 2)Modify it to fit the diffraction propagation formula
% us=us/Nx/psx/RI; vs=vs/Ny/psy/RI;
% HStack = exp(1i*lambda*pi*RI*(us.^2+vs.^2)* z);

% 2.1£©phase distribution under paraxial approximation (Equivalent to 1.2)
k0 = 2*pi/lambda;
[kx, ky] = meshgrid(cx, cy);
circle = (kx).^2+(ky).^2; 
R=Ny/2;
k = sqrt(kx.^2+ky.^2);
dkxy = NA/RI/floor(Nx/2);
sinthetaObj = k.*dkxy;
% sinthetaObj(sinthetaObj>1) = 1;
sinthetaObj(find(circle>R^2)) = 1;
phid = (sqrt(-1)*k0*RI).*(1-(sinthetaObj).^2/2); 
% phid = (sqrt(-1)*k0*RI).*((sinthetaObj).^2/2);
HStack = exp(phid * z); 

% %2.2£©cosine phase profile for spherical aberration compensation: fai=kz*cos(theta)
% costhetaObj = eps+sqrt(1-(sinthetaObj.^2));
% phid = (sqrt(-1)*k0*RI).*costhetaObj;
% HStack = exp(phid * z);

end
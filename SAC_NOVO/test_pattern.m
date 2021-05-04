%%% main program
%%% JCC 05/04/21
%%% In the subprogram, the one with "cos" in the name means program with
%%% spherical aberration correction, and the one without "cos" in the name means program without spherical aberration correction.
%% parameter definition
System.verbose=1;           % 1 or 0    Set this value to 1 to display activity, 0 otherwise
System.lambda = 1.04e-6;    % wavelength (m)
System.psSLM = 20e-6;       % SLM pixel size (m)
System.Nx = 600;            % Number of pixels in X direction (int)
System.Ny = 600;            % Number of pixels in Y direction (int)
System.useGPU = 0;          % 1 or 0    Use GPU to accelerate computation
System.maxiter = 100;       % Number of iterations (int)    
System.GSoffset = 1e-6;     % float>0   Regularization constant to allow low light background in 3D Gerchberg Saxton algorithms
NormOptions.HighThreshold = 0.9999;
NormOptions.LowThreshold = 0.001;
System.intensity = 1;       % Normalized total energy
System.source = sqrt(System.intensity)*(1/(System.Nx* System.Ny))*ones(System.Nx, System.Ny); %Definition of incident light field
x1=-System.Nx/2:System.Nx/2-1;   % 2D coordinate system definition of SLM
y1=-System.Ny/2:System.Ny/2-1;   % 2D coordinate system definition of SLM
RR=floor(System.Nx/2);           % The modulation boundary determined by the number of SLM
[x,y]=meshgrid(x1,y1);
circle = (x).^2+(y).^2;          
System.source(find(circle>RR^2))=0;  
System.focal_L1 = 0.4;        % Focal length of the first lens after slm (m)
System.focal_L2 = 0.5;        % Focal length of the second lens after slm (m)
System.focal_obj = 0.007;     % Focal length of the objective (m)
System.ObjNA = 1.05;          % NA of the objective
System.ObjRI = 1.33;          % The refractive index of the objective
Fx = System.lambda * System.focal_obj * System.focal_L1/ System.focal_L2 / System.psSLM; % Lateral field of view in image space
Fz = System.lambda * (System.focal_obj * System.focal_L1 / System.focal_L2)^2 / System.psSLM /(System.psSLM*System.Nx);  % Axial field of view in image space
delta_x = Fx/System.Nx;    
delta_dz = 1e-6;               % axial interval (m)

%% definition of target intensity distribution
% Here, Fig.5 is taken as an example: 
% bright patterns: radius: 5um, tilt_x:[30um 10um 30um], tilt_y:0, depth:[-50um 0 50um]
% dark patterns: radius: 5um, tilt_x:[20um 20um], tilt_y:0, depth:[-25um 25um]
R_t = 5e-6;                   % radius of target pattern(s) (m)
% Depths = 0e-6;                % depth of target pattern(s) (m)
Depths = -50e-6:25e-6:50e-6;  % depth of target pattern(s) (m)
tilt_x = [30e-6 20e-6 10e-6 20e-6 30e-6];     
tilt_y = 0e-6;       
TargetIntensity_bright = zeros(System.Nx,System.Ny,length(Depths)); 
TargetIntensity_dark = zeros(System.Nx,System.Ny,length(Depths));

for n = 1:length(Depths)    
    tiltX = round(tilt_x(n)/delta_x);
    tiltY = round(tilt_y/delta_x);
    circle2 = (x-tiltX).^2+(y-tiltY).^2;
    r_num = round(R_t/delta_x);  
    TargetIntensity_tem = zeros(System.Nx,System.Ny);
    TargetIntensity_tem(find(circle2<=r_num^2))=1;
    if mod(n,2)==1
        TargetIntensity_bright(:,:,n) = TargetIntensity_tem; % bright targets definition
    else
        TargetIntensity_dark(:,:,n) = TargetIntensity_tem;  % dark targets definition 
    end
end

%% Transfer function generation at differnt depth
[HStacks_com] = function_Hstacks_cos(System,Depths); % transfer function with spherical aberration compensation 
[HStacks] = function_Hstacks(System,Depths); % transfer function without spherical aberration compensation

%% phase searching
Hologram_sacnovo= function_NOVO_CGH_binary(System,HStacks_com,TargetIntensity_bright,TargetIntensity_dark,Depths,NormOptions);  % SAC-NOVO_Binary
Hologram_novo= function_NOVO_CGH_binary(System,HStacks,TargetIntensity_bright,TargetIntensity_dark,Depths,NormOptions);  % NOVO_Binary
phase_sacnovo = Hologram_sacnovo.phase+pi;
phase_novo = Hologram_novo.phase+pi;

%% 3D intensity distribution generation
Depths2 = Depths(1)-50e-6:delta_dz:Depths(length(Depths))+50e-6; % The axial range (m)
[ HStacks_com2 ] = function_Hstacks_cos( System,Depths2 );% transfer function with spherical aberration compensation 
Intenity_stack_sacnovo = function_VolumeIntensity( System, phase_sacnovo, HStacks_com2 ); % single photon 3D intensity distribution SAC-NOVO
Intenity_stack_novo = function_VolumeIntensity( System, phase_novo, HStacks_com2 ); % single photon 3D intensity distribution NOVO
Intenity_stack_sacnovo2 = Intenity_stack_sacnovo.^2; % two photon 3D intensity distribution SAC-NOVO
Intenity_stack_novo2 = Intenity_stack_novo.^2; % two photon 3D intensity distribution NOVO

%% GS pattern generation
% target definition
R_t = 5e-6;                   % radius of target pattern(s) (m)
Depths = -50e-6:50e-6:50e-6;  % depth of target pattern(s) (m)
tilt_x = [30e-6 10e-6 30e-6];     
tilt_y = 0e-6;       
TargetIntensity = zeros(System.Nx,System.Ny,length(Depths)); 

for n = 1:length(Depths)    
    tiltX = round(tilt_x(n)/delta_x);
    tiltY = round(tilt_y/delta_x);
    circle2 = (x-tiltX).^2+(y-tiltY).^2;
    r_num = round(R_t/delta_x);  
    TargetIntensity_tem = zeros(System.Nx,System.Ny);
    TargetIntensity_tem(find(circle2<=r_num^2))=1;
    TargetIntensity(:,:,n) = TargetIntensity_tem; % bright targets definition
end

% Transfer function generation
[HStacks_com] = function_Hstacks_cos(System,Depths); % transfer function with spherical aberration compensation 

% phase searching
Hologram_gs= function_globalGS(System,HStacks_com,TargetIntensity); % GS
% Hologram_VarI = function_NOVO_CGH_VarI(System,HStacks_com,TargetIntensity,NormOptions,Depths); % SAC-NOVO_VarI
phase_gs = Hologram_gs.phase+pi;
% phase_VarI = Hologram_VarI.phase+pi;

% 3D intensity distribution generation
Intenity_stack_gs = function_VolumeIntensity( System, phase_gs, HStacks_com2 ); % single photon 3D intensity distribution GS
Intenity_stack_gs2 = Intenity_stack_gs.^2; % two photon 3D intensity distribution GS

%% figure display
% phase
phase_mask = ones(System.Nx, System.Ny);
phase_mask(find(circle>RR^2))=0;
figure,
subplot(1,3,1),imshow(mod(phase_gs,2*pi).*phase_mask,[min(phase_gs(:)),max(phase_gs(:))]),title('GS phase');
subplot(1,3,2),imshow(mod(phase_novo,2*pi).*phase_mask,[min(phase_novo(:)),max(phase_novo(:))]),title('NOVO phase');
subplot(1,3,3),imshow(mod(phase_sacnovo,2*pi).*phase_mask,[min(phase_sacnovo(:)),max(phase_sacnovo(:))]),title('SAC-NOVO phase');

% intensity at XZsection
I_gs_xz = rot90(squeeze(Intenity_stack_gs2(301,300:390,:)),3);
figure,subplot(1,3,1),imshow(I_gs_xz,[min(I_gs_xz(:)),0.5*max(I_gs_xz(:))]),title('GS Intensity');

I_novo_xz = rot90(squeeze(Intenity_stack_novo2(301,300:390,:)),3);
subplot(1,3,2),imshow(I_novo_xz,[min(I_novo_xz(:)),0.5*max(I_novo_xz(:))]),title('NOVO Intensity');
I_sacnovo_xz = rot90(squeeze(Intenity_stack_sacnovo2(301,300:390,:)),3);
subplot(1,3,3),imshow(I_sacnovo_xz,[min(I_sacnovo_xz(:)),0.5*max(I_sacnovo_xz(:))]),title('SAC-NOVO Intensity');



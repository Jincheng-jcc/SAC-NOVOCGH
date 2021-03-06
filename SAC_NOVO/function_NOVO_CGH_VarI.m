function [ NOVOCGH ] = function_NOVO_CGH_VarI( System, HStacks, Masks,NormOptions,depths)
% This computes NOVO-CGH holgorams with binary targets. 
if System.verbose == 1
    disp(['NOVO-CGH with Variable Intensity Targets + Threshold based cost function, computation begins']);
end;
%Initialize phase guess with a superposition
vv = System.verbose;
System.verbose = 0;
[ Superposition ] = function_Superposition( System,HStacks,Masks);


%Calculate intensity in reconstruction to normalize thresholds
[ Superposition.Intensity_Reconstruction ] = function_VolumeIntensity(System, Superposition.phase,HStacks );
System.verbose = vv;
pct9999 = prctile(Superposition.Intensity_Reconstruction(:), 99.99);
prctile(Superposition.Intensity_Reconstruction(:), 99.9);
pct100 = max(Superposition.Intensity_Reconstruction(:));
thresholdh = NormOptions.HighThreshold * (pct9999 + pct100)/2;
thresholdl = NormOptions.LowThreshold * thresholdh;

% Define functions for mask kernel, and optimization local paramenters
kernelfun = @(i1, i2) HStacks(:,:,i1:i2);
maskfun = @(i1, i2) Masks(:,:,i1:i2);

f = @(x)function_FunObj_VarI( x, System.source, depths, System.Nx, System.Ny, thresholdh, thresholdl, maskfun, kernelfun, System.useGPU);
if System.verbose == 1; tic;end;
if System.verbose==1 ; displayoption = 'iter'; else; displayoption = 'off'; end
matlab_options = optimoptions('fmincon','GradObj','on', 'display', displayoption, ...
    'algorithm','interior-point','Hessian','lbfgs', 'MaxFunEvals', 500, 'MaxIter', System.maxiter,...
    'TolX', 1e-30, 'TolFun', 1e-18 , 'OutputFcn', @myoutput);
% options = trainingOptions('adam', ...
%     'InitialLearnRate',3e-4, ...
%     'SquaredGradientDecayFactor',0.99, ...
%     'MaxEpochs',20, ...
%     'MiniBatchSize',64, ...
%     'Plots','training-progress');
lb = -inf(System.Nx*System.Ny, 1);
ub = inf(System.Nx*System.Ny, 1);

%Optimize 
history_phase = [];
nn = 0;
[phase,fval] = fmincon(f,Superposition.phase,[],[],[],[],lb,ub,[],matlab_options);
    function stop = myoutput(phase,optimvalues,state)
        stop = false;
        if isequal(state,'iter')
          nn = nn+1;
          history_phase(:,:,nn) = phase;
        end
    end
%Reformat result to extract hologram
phase = reshape(phase, [System.Nx, System.Ny]);
phase = mod(phase, 2*pi) - pi;
NOVOCGH.phase=phase;

if System.verbose == 1
    t = toc; 
    disp('NOVO-CGH with Variable Intensity Targets + Threshold based cost function')
    disp('And with Dark ROI enforcement')
    disp(['- Completed in ' int2str(t) ' seconds !']);
end;
end


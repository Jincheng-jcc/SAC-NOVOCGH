%%% series of target intensity (ims), dark voxels (darkims), and depths (imdepths) for the desired format. 
%%% 层深的单位：m
function [ims_all,dims,imdepths,sphere,ims_slice,sphere_mask_whole ] = function_extended_pattern(System)

% psXHolograph = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Nx / System.focal_tube;      % Pixel Size (resolution) at the scattered 3D region
% psYHolograph = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Ny / System.focal_tube;      % Pixel Size (resolution) at the scattered 3D region
% 
% imdepths = [-15e-6,0,15e-6];
% % imdepths = 0;
% dims = ones(System.Nx,System.Ny,length(imdepths)); 
% dims(:,:,2)=0;
% ims_all = zeros(System.Nx,System.Ny,length(imdepths)); %3D目标包含的体素
% 
% x1=-System.Nx/2:System.Nx/2-1;   %把圆心变到矩阵的中间
% y1=-System.Ny/2:System.Ny/2-1;
% [x,y]=meshgrid(x1,y1);
% tem_pattern = zeros(System.Nx,System.Ny);
% circle=(x-x1(System.Nx/2)).^2+(y-y1(System.Ny/2)).^2;   %计算出每一点到圆心的距离的平方
% sphere =[];
%% 1.单层扩展pattern生成
% pattern_size = 10e-6; %扩展pattern的大小
% r_num = round(pattern_size/psXHolograph);
% tem_pattern(find(circle<=(r_num^2)))=1;
% ims_all(:,:,2) = tem_pattern;
% % ims_all(:,:,1) = tem_pattern;

%% 2.生成环形pattern
% pattern_size1 = 10e-6; %外径尺寸
% pattern_size2 = 8e-6; %内径尺寸
% r_num1 = round(pattern_size1/psXHolograph);
% r_num2 = round(pattern_size2/psXHolograph);
% tem_pattern = ones(System.Nx,System.Ny);
% tem_pattern(find((circle>=(r_num1^2))))=0;
% tem_pattern(find((circle<=(r_num2^2))))=0;
% ims_all(:,:,2) = tem_pattern;
% % ims_all(:,:,1) = tem_pattern;
% ims = ims_all;

%% 3.3D multi-sphere pattern
% % 定义3D坐标系
% Fx = System.lambda  * System.focal_obj * System.focal_SLM / System.focal_tube / System.psSLM /2;
% Fz = System.lambda  * (System.focal_obj * System.focal_SLM / System.focal_tube)^2 / System.psSLM /(System.psSLM*System.Nx)/2;
% delta_z = 1e-6;
% delta_x = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Nx / System.focal_tube; 
% % Nz = round(2*Fz/delta_z);
% Nz = round(60e-6/delta_z);
% x1=-System.Nx/2:System.Nx/2-1;   %把圆心变到矩阵的中间
% y1=-System.Ny/2:System.Ny/2-1;
% z1=-Nz/2:Nz/2-1;
% [x,y]=meshgrid(x1,y1);
% dims=zeros(System.Nx,System.Ny,Nz); 
% imdepths=delta_z*z1;
% % 定义pattern大小&数目&位置
% sphere.R = 5e-6; %球体半径
% % sphere.NUM = 15; %球体数目
% sphere.NUM = 5; %球体数目
% sphere.center = zeros(sphere.NUM,3); %球心位置
% sphere.delta_z = delta_z;
% n=1;
% r_num_z = round(sphere.R/delta_z);
% r_num_x = round(sphere.R/delta_x);
% % 缩小横向定位范围
% r_num_x = round(10*sphere.R/delta_x);
% % while n<=sphere.NUM
%     % 随机pattern中心
%     sphere.center(n,1)=randi([r_num_x,System.Nx-r_num_x],1,1);
%     sphere.center(n,2)=randi([r_num_x,System.Ny-r_num_x],1,1);
%     sphere.center(n,3)=randi([r_num_z,Nz-r_num_z],1,1);
% 
% % % 单pattern中心定位
% %     sphere.center(n,1)=round(System.Nx/2);
% %     sphere.center(n,2)=round(System.Ny/2);
% %     sphere.center(n,3)=round(Nz/2);
% %     flag = 0;
% %     for i = 1:n-1
% %         d_x = abs(sphere.center(n,1)-sphere.center(i,1));
% %         d_y = abs(sphere.center(n,2)-sphere.center(i,2));
% %         d_z = abs(sphere.center(n,3)-sphere.center(i,3));
% % %         if (d_x<2*r_num_x)||(d_y<2*r_num_x)||(d_z<2*r_num_z)
% %         if (d_x<2*round(sphere.R/delta_x))||(d_y<2*round(sphere.R/delta_x))
% %             flag=1;
% %             nn = n
% %             break
% %         end
% %         
% %     end
% %     if flag==0
% %         n=n+1;
% %         n
% %     end
% % end
% % 定中心坐标球形 (轴向)
% sphere.center(:,1)= [round(System.Nx/4),round(3*System.Nx/4),round(System.Nx/2),round(System.Nx/4),round(3*System.Nx/4)];
% sphere.center(:,2)= [round(System.Ny/4),round(System.Ny/4),round(System.Ny/2),round(3*System.Ny/4),round(3*System.Ny/4)];
% sphere.center(:,3)= [round(Nz/6),round(Nz/3),round(Nz/2),round(2*Nz/3),round(5*Nz/6)];
% % % 定中心坐标球形 (横向)
% % sphere.center(:,1)= [round(System.Nx/4),round(3*System.Nx/4),round(System.Nx/2),round(System.Nx/4),round(3*System.Nx/4)];
% % sphere.center(:,2)= [round(System.Ny/4),round(System.Ny/4),round(System.Ny/2),round(3*System.Ny/4),round(3*System.Ny/4)];
% % sphere.center(:,3)= [round(Nz/2),round(Nz/2),round(Nz/2),round(Nz/2),round(Nz/2)];
% % 生成对应的3Dpattern
% sphere_mask = zeros(System.Nx,System.Ny,Nz); 
% for k = 1:sphere.NUM
% sphere_mask_tem = zeros(System.Nx,System.Ny,Nz); 
% zh1 = r_num_z+sphere.center(k,3);
% zh2 = sphere.center(k,3)-r_num_z; % 轴向采样范围下限
% % zh1 = sphere.center(k,3);
% % zh2 = sphere.center(k,3); % 生成圆盘目标函数
% circle = (x-x1(sphere.center(k,1))).^2+(y-y1(sphere.center(k,2))).^2;   %计算出每一点到圆心的距离的平方
% 
% for zh = zh2:zh1
%     r = sqrt(sphere.R^2-(abs(zh-sphere.center(k,3))*delta_z)^2);
%     r_num = round(r/delta_x);
%     circ_mask = zeros(System.Nx,System.Ny);
%     circ_mask(find(circle<r_num^2))=1;  %找到圆内的元素，赋值为1
%     sphere_mask_tem(:,:,zh) = circ_mask;    
% end
% 
% sphere_mask = sphere_mask|sphere_mask_tem;
% sphere_mask = double(sphere_mask);
% 
% end
% 
% ims_all = sphere_mask;
% sphere_mask_index = uint8(sphere_mask*255);
% % saveastiff(uint8(sphere_mask_index), 'sphere_mask_15_beads_D10um.tif');
% saveastiff(uint8(sphere_mask_index), 'define_mask_5beads_D10um_30um.tif');
%% 4.密集pattern生成
% % 定义3D坐标系
% Fx = System.lambda  * System.focal_obj * System.focal_SLM / System.focal_tube / System.psSLM /2;
% Fz = System.lambda  * (System.focal_obj * System.focal_SLM / System.focal_tube)^2 / System.psSLM /(System.psSLM*System.Nx)/2;
% delta_z = 1e-6;
% delta_x = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Nx / System.focal_tube; 
% Nz = round(2*Fz/delta_z);
% x1=-System.Nx/2:System.Nx/2-1;   %把圆心变到矩阵的中间
% y1=-System.Ny/2:System.Ny/2-1;
% z1=-Nz/2:Nz/2-1;
% [x,y]=meshgrid(x1,y1);
% dims=zeros(System.Nx,System.Ny,Nz); 
% imdepths=delta_z*z1;
% % 定义pattern大小&数目&位置
% sphere.R = 5e-6; %球体半径
% % sphere.NUM = 15; %球体数目
% sphere.NUM = 8; %球体数目
% sphere.center = zeros(sphere.NUM,3); %球心位置
% sphere.delta_z = delta_z;
% distance_x = 1.05* sphere.R;
% distance_y = 1.05* sphere.R*sqrt(3);
% r_num_z = round(sphere.R/delta_z);
% r_num_x = round(distance_x/delta_x);
% r_num_y = round(distance_y/delta_x);
% % 定义密集排布Pattern [参考Eirini Papagiakoumou et al. 2018Two-Photon Optogenetics by Computer-Generated Holography]
% sphere.center(:,1)= [round(System.Nx/2)-4*r_num_x,round(System.Nx/2),round(System.Nx/2)+r_num_x,round(System.Nx/2)+r_num_x,round(System.Nx/2)+2*r_num_x,round(System.Nx/2)+3*r_num_x,round(System.Nx/2)+3*r_num_x,round(System.Nx/2)+4*r_num_x];
% sphere.center(:,2)= [round(System.Ny/2),round(System.Ny/2),round(System.Ny/2)+r_num_y,round(System.Ny/2)-r_num_y,round(System.Ny/2),round(System.Ny/2)+r_num_y,round(System.Ny/2)-r_num_y,round(System.Ny/2)];
% sphere.center(:,3)=ones(sphere.NUM,1)*round(Nz/2);
% 
% 
%     sphere.center(n,1)=randi([r_num_x,System.Nx-r_num_x],1,1);
%     sphere.center(n,2)=randi([r_num_x,System.Ny-r_num_x],1,1);
%     sphere.center(n,3)=randi([r_num_z,Nz-r_num_z],1,1);
% 
% % 生成对应的3Dpattern
% sphere_mask = zeros(System.Nx,System.Ny,Nz); 
% for k = 1:sphere.NUM
% sphere_mask_tem = zeros(System.Nx,System.Ny,Nz); 
% zh1 = r_num_z+sphere.center(k,3);
% zh2 = sphere.center(k,3)-r_num_z; % 轴向采样范围下限
% circle = (x-x1(sphere.center(k,1))).^2+(y-y1(sphere.center(k,2))).^2;   %计算出每一点到圆心的距离的平方
% 
% for zh = zh2:zh1
%     r = sqrt(sphere.R^2-(abs(zh-sphere.center(k,3))*delta_z)^2);
%     r_num = round(r/delta_x);
%     circ_mask = zeros(System.Nx,System.Ny);
%     circ_mask(find(circle<=r_num^2))=1;  %找到圆内的元素，赋值为1
%     sphere_mask_tem(:,:,zh) = circ_mask;    
% end
% 
% sphere_mask = sphere_mask|sphere_mask_tem;
% sphere_mask = double(sphere_mask);
% 
% end
% ims_all = sphere_mask;
% sphere_mask_index = uint8(sphere_mask*255);
% % saveastiff(uint8(sphere_mask_index), 'sphere_mask_15_beads_D10um.tif');
% saveastiff(uint8(sphere_mask_index), 'sphere_mask_dense_beads_D10um.tif');

%% 5.单层pattern生成
% delta_x = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Nx / System.focal_tube; 
% % Nz = round(2*Fz/delta_z);
% x1=-System.Nx/2:System.Nx/2-1;   %把圆心变到矩阵的中间
% y1=-System.Ny/2:System.Ny/2-1;
% [x,y]=meshgrid(x1,y1);
% % imdepths=[-40e-6 -20e-6 0 20e-6 40e-6];
% % imdepths=[0 0 0 0 0];
% delta_z = System.lambda * (System.focal_obj * System.focal_SLM / System.focal_tube / System.psSLM  / System.Nx)^2; %由SLM决定的轴向分辨率
% Fz = System.lambda * (System.focal_obj * System.focal_SLM / System.focal_tube)^2 / System.psSLM /(System.psSLM*System.Nx)/2;  %由SLM决定的轴向视场
% Nz=round(2*Fz/delta_z);
% z1=-Nz/2:Nz/2-1;
% z1 = z1*delta_z;
% % imdepths=z1([round(Nz/6),round(Nz/3),round(Nz/2),round(2*Nz/3),round(5*Nz/6)]); %均匀轴向间隔
% % imdepths=z1([100 200 300 400 500]); %非均匀轴向间隔
% % imdepths=z1([260 300 340]);
% imdepths=z1([300]);
% nz = length(imdepths);
% dims=zeros(System.Nx,System.Ny,nz); 
% % 定义pattern大小&数目&位置
% sphere.R = 5e-6; %球体半径
% % R = [30e-6,10e-6,30e-6];
% R = [10e-6];
% % sphere.NUM = 15; %球体数目
% sphere.NUM = 1; %球体数目
% sphere.center = zeros(sphere.NUM,3); %球心位置
% 
% % % 定中心坐标球形 
% % sphere.center(:,1)= [round(System.Nx/4),round(3*System.Nx/4),round(System.Nx/2),round(System.Nx/4),round(3*System.Nx/4)];
% % sphere.center(:,2)= [round(System.Ny/4),round(System.Ny/4),round(System.Ny/2),round(3*System.Ny/4),round(3*System.Ny/4)];
% % sphere.center(:,3)= [round(Nz/2),round(Nz/2),round(Nz/2),round(Nz/2),round(Nz/2)]; (横向)
% % sphere.center(:,1)= [round(System.Nx*0.3),round(0.35*System.Nx),round(0.6*System.Nx),round(0.75*System.Nx),round(0.8*System.Nx)];
% % sphere.center(:,2)= [round(0.2*System.Ny),round(0.8*System.Ny),round(0.7*System.Ny),round(0.4*System.Ny),round(0.6*System.Ny)];
% % sphere.center(:,1)= [round(System.Nx*0.2),round(0.35*System.Nx),round(0.5*System.Nx),round(0.65*System.Nx),round(0.8*System.Nx)];
% % sphere.center(:,2)= [round(0.2*System.Ny),round(0.35*System.Ny),round(0.5*System.Ny),round(0.65*System.Ny),round(0.8*System.Ny)];
% % sphere.center(:,3)= [1 2 3 4 5]; %直接给出球心轴向坐标
% sphere.center(:,1)= [round(System.Nx*0.5)];
% sphere.center(:,2)= [round(0.5*System.Ny)];
% sphere.center(:,3)= [1]; %直接给出球心轴向坐标
% % sphere.center(:,3)= [round(Nz/6),round(Nz/3),round(Nz/2),round(2*Nz/3),round(5*Nz/6)];  %（轴向）
% % 生成对应的3Dpattern
% sphere_mask = zeros(System.Nx,System.Ny,nz); 
% sphere_mask_slice = zeros(System.Nx,System.Ny); 
% % weight_intensity = ones(sphere.NUM,1);
% for zh = 1:length(imdepths)
% k = 1;
% circle = (x-x1(sphere.center(k,1))).^2+(y-y1(sphere.center(k,2))).^2;   %计算出每一点到圆心的距离的平方
% %     r = sqrt(sphere.R^2-(imdepths(zh))^2);
% %       r = sphere.R;
%     r = R(zh);
% %     r = delta_x;
%     r_num = round(r/delta_x);
%     circ_mask = zeros(System.Nx,System.Ny);
%     circ_mask(find(circle<=r_num^2))=1;  %找到圆内的元素，赋值为1
% %     circ_mask = circ_mask.*weight_intensity(zh);
% %     sphere_mask(:,:,zh) = sphere_mask(:,:,zh)|circ_mask;
%     sphere_mask(:,:,zh) = circ_mask;
%     sphere_mask_slice =sphere_mask_slice|sphere_mask(:,:,zh);
% end
% sphere_mask = double(sphere_mask);
% ims_all = sphere_mask;
% % ims_slice = double(sphere_mask_slice);
% ims_slice = ims_all;
% sphere_mask_index = uint8(sphere_mask*255);
% % saveastiff(uint8(sphere_mask_index), 'sphere_mask_15_beads_D10um.tif');
% % saveastiff(uint8(sphere_mask_index),['1beads_D',num2str(2*sphere.R*1e6),'um_depths',num2str(imdepths*1e6),'um.tif']);
% saveastiff(uint8(sphere_mask_index),['1beads_D',num2str(2*sphere.R*1e6),'um_depths_-10 0 10um.tif']);

%% 6. 3*3*3 晶格结构扩展pattern
% delta_x = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Nx / System.focal_tube; 
% x1=-System.Nx/2:System.Nx/2-1;   %把圆心变到矩阵的中间
% y1=-System.Ny/2:System.Ny/2-1;
% [x,y]=meshgrid(x1,y1);
% delta_z = System.lambda * (System.focal_obj * System.focal_SLM / System.focal_tube / System.psSLM  / System.Nx)^2; %由SLM决定的轴向分辨率
% Fz = System.lambda * (System.focal_obj * System.focal_SLM / System.focal_tube)^2 / System.psSLM /(System.psSLM*System.Nx)/2;  %由SLM决定的轴向视场
% Nz=round(2*Fz/delta_z);
% z1=-Nz/2:Nz/2-1;
% z1 = z1*delta_z;
% % zint = [1 300 600];
% % zint = [100 300 450];
% zint = [300];
% imdepths=z1(zint); %非均匀轴向间隔
% nz = length(imdepths);
% dims=zeros(System.Nx,System.Ny,nz); 
% % 定义pattern大小&数目&位置
% sphere.R = 5e-6; %球体半径
% sphere.Rint = 3e-6; % 环形球内径
% % sphere.NUM = 15; %球体数目
% sphere.NUM_perslice=9;
% % sphere.NUM_perslice=1;
% sphere.NUM = sphere.NUM_perslice*nz; %球体数目
% sphere.center = zeros(sphere.NUM,3); %球心位置
% xs= [round(System.Nx*0.33),round(0.5*System.Nx),round(0.67*System.Nx)];
% ys= [round(System.Ny*0.33),round(0.5*System.Ny),round(0.67*System.Ny)];
% % xs= [round(0.5*System.Nx)];
% % ys= [round(0.5*System.Ny)];
% [Xs,Ys]=meshgrid(xs,ys);
% Xs = reshape(Xs,[1,sphere.NUM_perslice]);
% Ys = reshape(Ys,[1,sphere.NUM_perslice]);
% sphere.center(:,1) = repmat(Xs,1,nz);
% sphere.center(:,2)=  repmat(Ys,1,nz);
% % sphere.center(:,3)= reshape(([1 2 3]'*ones(1,9))',[1,sphere.NUM]); %直接给出球心轴向坐标
% sphere.center(:,3)= reshape(([1]'*ones(1,9))',[1,sphere.NUM]); %直接给出球心轴向坐标
% % sphere.center(:,3)= reshape(([1]'*ones(1))',[1,sphere.NUM]); %直接给出球心轴向坐标
% sphere.mask = zeros(System.Nx,System.Ny,nz);
% % sphere.center(:,3)= [round(Nz/6),round(Nz/3),round(Nz/2),round(2*Nz/3),round(5*Nz/6)];  %（轴向）
% % 生成对应的3Dpattern
% sphere_mask = zeros(System.Nx,System.Ny,nz); 
% sphere_mask_slice = zeros(System.Nx,System.Ny); 
% % weight_intensity = ones(sphere.NUM,1);
% for zh = 1:length(imdepths)
%     for k = 1:sphere.NUM_perslice
%       kk = (zh-1)*sphere.NUM_perslice+k;
%       circle = (x-x1(sphere.center(kk,1))).^2+(y-y1(sphere.center(kk,2))).^2;   %计算出每一点到圆心的距离的平方
%       r_num = round(sphere.R/delta_x);
%       r_num2 = round(sphere.Rint/delta_x);
%       circ_mask = zeros(System.Nx,System.Ny);
%       circ_mask2 = zeros(System.Nx,System.Ny);
%       circ_mask(find(circle<=r_num^2))=1;  %找到圆内的元素，赋值为1
%       %环形圆盘
%       circ_mask2(find(circle<=r_num2^2))=1;
%       circ_mask=circ_mask-circ_mask2;  %找到环内的元素，赋值为1
%       sphere_mask(:,:,zh) = sphere_mask(:,:,zh)|circ_mask;
%       sphere_mask_slice =sphere_mask_slice|sphere_mask(:,:,zh);
%       sphere.mask(:,:,kk) = circ_mask;
%     end
% end
% sphere_mask = double(sphere_mask);
% ims_all = sphere_mask;
% % ims_slice = double(sphere_mask_slice);
% ims_slice = ims_all;
% sphere_mask_index = uint8(sphere_mask*255);
% saveastiff(uint8(sphere_mask_index), 'sphere_mask_27_beads_D10um.tif');
% saveastiff(uint8(sphere_mask_index), '1beads_D60um_1depth_1slice.tif');

% % 扩展pattern完整版
% Nz=round(2*Fz/delta_z)+100;
% zint = [1 300 600]+50;
% r_num_z = round(sphere.R/delta_z);
% sphere_mask_whole = zeros(System.Nx,System.Ny,Nz); 
% for k = 1:sphere.NUM
% zh1 = r_num_z+zint(sphere.center(k,3));
% zh2 = zint(sphere.center(k,3))-r_num_z; % 轴向采样范围下限
% circle = (x-x1(sphere.center(k,1))).^2+(y-y1(sphere.center(k,2))).^2;   %计算出每一点到圆心的距离的平方
% for zh = zh2:zh1
%     r = sqrt(sphere.R^2-(abs(zh-zint(sphere.center(k,3)))*delta_z)^2);
%     r_num = round(r/delta_x);
%     circ_mask = zeros(System.Nx,System.Ny);
%     circ_mask(find(circle<=r_num^2))=1;  %找到圆内的元素，赋值为1
%     sphere_mask_whole(:,:,zh) = sphere_mask_whole(:,:,zh)|circ_mask;
% end
% end
% sphere_mask_whole = double(sphere_mask_whole);
% ims_all = sphere_mask_whole;
% sphere_mask_index = uint8(sphere_mask_whole*255);
% % saveastiff(uint8(sphere_mask_index), 'sphere_mask_15_beads_D10um.tif');
% saveastiff(uint8(sphere_mask_index), '27beads_whole_D10um_3depths_3slice.tif');
%% 7.扩展pattern轴向压窄
% delta_x = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Nx / System.focal_tube; 
% % Nz = round(2*Fz/delta_z);
% x1=-System.Nx/2:System.Nx/2-1;   %把圆心变到矩阵的中间
% y1=-System.Ny/2:System.Ny/2-1;
% [x,y]=meshgrid(x1,y1);
% % imdepths=[-40e-6 -20e-6 0 20e-6 40e-6];
% % imdepths=[0 0 0 0 0];
% delta_z = System.lambda * (System.focal_obj * System.focal_SLM / System.focal_tube / System.psSLM  / System.Nx)^2; %由SLM决定的轴向分辨率
% Fz = System.lambda * (System.focal_obj * System.focal_SLM / System.focal_tube)^2 / System.psSLM /(System.psSLM*System.Nx)/2;  %由SLM决定的轴向视场
% Nz=round(2*Fz/delta_z);
% z1=-Nz/2:Nz/2-1;
% z1 = z1*delta_z;
% % imdepths=z1([round(Nz/6),round(Nz/3),round(Nz/2),round(2*Nz/3),round(5*Nz/6)]); %均匀轴向间隔
% % imdepths=z1([100 200 300 400 500]); %非均匀轴向间隔
% % imdepths=z1([280 300 320]);
% imdepths=z1([10:20:280,300,320:20:590]);
% % imdepths=z1([300]);
% nz = length(imdepths);
% dims=zeros(System.Nx,System.Ny,nz); 
% % 定义pattern大小&数目&位置
% sphere.R = 5e-6; %球体半径
% % R = [30e-6,10e-6,30e-6];
% R = [10e-6];
% % sphere.NUM = 15; %球体数目
% sphere.NUM = 1; %球体数目
% sphere.center = zeros(sphere.NUM,3); %球心位置
% 
% % % 定中心坐标球形 
% % sphere.center(:,1)= [round(System.Nx/4),round(3*System.Nx/4),round(System.Nx/2),round(System.Nx/4),round(3*System.Nx/4)];
% % sphere.center(:,2)= [round(System.Ny/4),round(System.Ny/4),round(System.Ny/2),round(3*System.Ny/4),round(3*System.Ny/4)];
% % sphere.center(:,3)= [round(Nz/2),round(Nz/2),round(Nz/2),round(Nz/2),round(Nz/2)]; (横向)
% % sphere.center(:,1)= [round(System.Nx*0.3),round(0.35*System.Nx),round(0.6*System.Nx),round(0.75*System.Nx),round(0.8*System.Nx)];
% % sphere.center(:,2)= [round(0.2*System.Ny),round(0.8*System.Ny),round(0.7*System.Ny),round(0.4*System.Ny),round(0.6*System.Ny)];
% % sphere.center(:,1)= [round(System.Nx*0.2),round(0.35*System.Nx),round(0.5*System.Nx),round(0.65*System.Nx),round(0.8*System.Nx)];
% % sphere.center(:,2)= [round(0.2*System.Ny),round(0.35*System.Ny),round(0.5*System.Ny),round(0.65*System.Ny),round(0.8*System.Ny)];
% % sphere.center(:,3)= [1 2 3 4 5]; %直接给出球心轴向坐标
% sphere.center(:,1)= [round(System.Nx*0.5)];
% sphere.center(:,2)= [round(0.5*System.Ny)];
% sphere.center(:,3)= [1]; %直接给出球心轴向坐标
% % sphere.center(:,3)= [round(Nz/6),round(Nz/3),round(Nz/2),round(2*Nz/3),round(5*Nz/6)];  %（轴向）
% % 生成对应的3Dpattern
% sphere_mask = zeros(System.Nx,System.Ny,nz); 
% sphere_mask_slice = zeros(System.Nx,System.Ny); 
% k = 1;
% circle = (x-x1(sphere.center(k,1))).^2+(y-y1(sphere.center(k,2))).^2;   %计算出每一点到圆心的距离的平方
% %     r = sqrt(sphere.R^2-(imdepths(zh))^2);
%       r = sphere.R;
% %     r = R(zh);
% %     r = delta_x;
%     r_num = round(r/delta_x);
%     circ_mask = zeros(System.Nx,System.Ny);
%     circ_mask(find(circle<=r_num^2))=1;  %找到圆内的元素，赋值为1
%     
% I_sum = sum(circ_mask(:));
% for m = 1:length(imdepths)
% sphere_mask(round(System.Nx/4):round(3:System.Nx/4),round(System.Ny/4):round(3:System.Ny/4),m) = ones(round(System.Nx/4):round(3:System.Nx/4),round(System.Ny/4):round(3:System.Ny/4)).*I_sum/(System.Nx*System.Ny/4);
% end
% sphere_mask(:,:,round(length(imdepths)/2)) = circ_mask;
% %     circ_mask = circ_mask.*weight_intensity(zh);
% %     sphere_mask(:,:,zh) = sphere_mask(:,:,zh)|circ_mask;
% % sphere_mask(:,:,1) = ones(System.Nx,System.Ny);
% % sphere_mask(:,:,3) = ones(System.Nx,System.Ny);
% 
% sphere_mask = double(sphere_mask);
% ims_all = sphere_mask;
% % ims_slice = double(sphere_mask_slice);
% ims_slice = ims_all;
% sphere_mask_index = uint8(sphere_mask*255);
% % saveastiff(uint8(sphere_mask_index), 'sphere_mask_15_beads_D10um.tif');
% % saveastiff(uint8(sphere_mask_index),['1beads_D',num2str(2*sphere.R*1e6),'um_depths',num2str(imdepths*1e6),'um.tif']);
% saveastiff(uint8(sphere_mask_index),['energylimit_1beads_D',num2str(2*sphere.R*1e6),'um_slice',num2str(imdepths),'.tif']);
%% 8. 3D环形pattern生成，与3D环形disk对比
delta_x = System.lambda * System.focal_SLM * System.focal_obj / System.psSLM  / System.Nx / System.focal_tube; 
% Nz = round(2*Fz/delta_z);
x1=-System.Nx/2:System.Nx/2-1;   %把圆心变到矩阵的中间
y1=-System.Ny/2:System.Ny/2-1;
[x,y]=meshgrid(x1,y1);
% imdepths=[-40e-6 -20e-6 0 20e-6 40e-6];
% imdepths=[0 0 0 0 0];
delta_z = System.lambda * (System.focal_obj * System.focal_SLM / System.focal_tube / System.psSLM  / System.Nx)^2; %由SLM决定的轴向分辨率
Fz = System.lambda * (System.focal_obj * System.focal_SLM / System.focal_tube)^2 / System.psSLM /(System.psSLM*System.Nx)/2;  %由SLM决定的轴向视场
Nz=round(2*Fz/delta_z);
z1=-Nz/2:Nz/2-1;
z1 = z1*delta_z;
% imdepths=z1(100:100:500);
% imdepths=z1([300]);R = [10e-6];
% sphere.NUM = 15; %球体数目
sphere.NUM = 5; %球体数目
imdepths=z1(ones(1,sphere.NUM)*300);
nz = length(imdepths);
dims=zeros(System.Nx,System.Ny,nz); 
% 定义pattern大小&数目&位置
sphere.R = 10e-6; %球体半径
sphere.R_in = 9e-6; %球体半径
% r_num = round(sphere.R/delta_x)+2;
% R = [30e-6,10e-6,30e-6];

sphere.center = zeros(sphere.NUM,3); %球心位置
% sphere.center(:,3)= [1 2 3 4 5]; %直接给出球心轴向坐标
% sphere.center(:,1)= [round(System.Nx/2)-4*r_num,round(System.Nx/2),round(System.Nx/2)+r_num,round(System.Nx/2)+r_num,round(System.Nx/2)+2*r_num,round(System.Nx/2)+3*r_num,round(System.Nx/2)+3*r_num,round(System.Nx/2)+4*r_num];
% sphere.center(:,1)= [round(System.Nx/2)-4*r_num,round(System.Nx/2),round(System.Nx/2)+r_num,round(System.Nx/2)+r_num,round(System.Nx/2)+2*r_num,round(System.Nx/2)+3*r_num,round(System.Nx/2)+3*r_num,round(System.Nx/2)+4*r_num];
% sphere.center(:,2)= [round(System.Ny/2),round(System.Ny/2),round(System.Ny/2)+2*r_num,round(System.Ny/2)-2*r_num,round(System.Ny/2),round(System.Ny/2)+2*r_num,round(System.Ny/2)-2*r_num,round(System.Ny/2)];

sphere.center(:,1)= [round(System.Nx/4),round(3*System.Nx/4),round(System.Nx/2),round(System.Nx/4),round(3*System.Nx/4)];
sphere.center(:,2)= [round(System.Ny/4),round(System.Ny/4),round(System.Ny/2),round(3*System.Ny/4),round(3*System.Ny/4)];
sphere.center(:,3)= 1:5; %直接给出球心轴向坐标
% 生成对应的3Dpattern
sphere_mask = zeros(System.Nx,System.Ny,nz); 
ims_slice = zeros(System.Nx,System.Ny); 
r_num = round(sphere.R/delta_x);
for k = 1:nz
    circle = (x-x1(sphere.center(k,1))).^2+(y-y1(sphere.center(k,2))).^2;   %计算出每一点到圆心的距离的平方
%     r = sphere.R;
%     r_num = round(r/delta_x);
    r_num2=round(sphere.R_in/delta_x);
    circ_mask = zeros(System.Nx,System.Ny);
    circ_mask(find(circle<=r_num^2))=1;  %找到圆内的元素，赋值为1
    circ_mask(find(circle<=r_num2^2))=0;  
    sphere_mask(:,:,k) = circ_mask;
    ims_slice = ims_slice|sphere_mask(:,:,k);
end
ims_slice= double(ims_slice);
sphere_mask = double(sphere_mask);
ims_all = sphere_mask;
% ims_slice = double(sphere_mask_slice);

sphere_mask_index = uint8(sphere_mask*255);
% saveastiff(uint8(sphere_mask_index), 'sphere_mask_15_beads_D10um.tif');
% saveastiff(uint8(sphere_mask_index),['1beads_D',num2str(2*sphere.R*1e6),'um_depths',num2str(imdepths*1e6),'um.tif']);
saveastiff(uint8(sphere_mask_index),['5dense_annular1um_D',num2str(2*sphere.R*1e6),'um_slice',num2str(length(imdepths)),'.tif']);
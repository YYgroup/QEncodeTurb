addpath('../plotTools/vtkwrite');
load('extractedData_3Dhomo.mat')

L2_max = max([lambda2_t200(:); lambda2_t800(:); lambda2_t1400(:); lambda2_t2000(:)]);

mkdir ../../results/VTKfiles/;

%% DNS
%DNS,t=200
L2 = lambda2_t200;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2_t200.vtk','structured_points','t200',L2);

%DNS,t=800
L2 = lambda2_t800;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2_t800.vtk','structured_points','t800',L2);

%DNS,t=1400
L2 = lambda2_t1400;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2_t1400.vtk','structured_points','t1400',L2);

%DNS,t=2000
L2 = lambda2_t2000;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2_t2000.vtk','structured_points','t2000',L2);

%% Chi=192
%C192,t=200
L2 = lambda2C192_t200;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C192_t200.vtk','structured_points','t200',L2);

%C192,t=800
L2 = lambda2C192_t800;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C192_t800.vtk','structured_points','t800',L2);

%C192,t=1400
L2 = lambda2C192_t1400;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C192_t1400.vtk','structured_points','t1400',L2);

%C192,t=2000
L2 = lambda2C192_t2000;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C192_t2000.vtk','structured_points','t2000',L2);

%% Chi=128
%C128,t=200
L2 = lambda2C128_t200;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C128_t200.vtk','structured_points','t200',L2);

%C128,t=800
L2 = lambda2C128_t800;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C128_t800.vtk','structured_points','t800',L2);

%C128,t=1400
L2 = lambda2C128_t1400;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C128_t1400.vtk','structured_points','t1400',L2);

%C128,t=2000
L2 = lambda2C128_t2000;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C128_t2000.vtk','structured_points','t2000',L2);

%% Chi=96
%C96,t=200
L2 = lambda2C96_t200;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C96_t200.vtk','structured_points','t200',L2);

%C96,t=800
L2 = lambda2C96_t800;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C96_t800.vtk','structured_points','t800',L2);

%C96,t=1400
L2 = lambda2C96_t1400;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C96_t1400.vtk','structured_points','t1400',L2);

%C96,t=2000
L2 = lambda2C96_t2000;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2C96_t2000.vtk','structured_points','t2000',L2);

%% N=88
%N88,t=200
L2 = lambda2N88_t200;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N88_t200.vtk','structured_points','t200',L2);

%N88,t=800
L2 = lambda2N88_t800;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N88_t800.vtk','structured_points','t800',L2);

%N88,t=1400
L2 = lambda2N88_t1400;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N88_t1400.vtk','structured_points','t1400',L2);

%N88,t=2000
L2 = lambda2N88_t2000;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N88_t2000.vtk','structured_points','t2000',L2);

%% N=70
%N70,t=200
L2 = lambda2N70_t200;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N70_t200.vtk','structured_points','t200',L2);

%N70,t=800
L2 = lambda2N70_t800;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N70_t800.vtk','structured_points','t800',L2);

%N70,t=1400
L2 = lambda2N70_t1400;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N70_t1400.vtk','structured_points','t1400',L2);

%N70,t=2000
L2 = lambda2N70_t2000;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N70_t2000.vtk','structured_points','t2000',L2);

%% N=60
%N60,t=200
L2 = lambda2N60_t200;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N60_t200.vtk','structured_points','t200',L2);

%N60,t=800
L2 = lambda2N60_t800;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N60_t800.vtk','structured_points','t800',L2);

%N60,t=1400
L2 = lambda2N60_t1400;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../results/VTKfiles/lambda2N60_t1400.vtk','structured_points','t1400',L2);

%N60,t=2000
L2 = lambda2N60_t2000;
L2(L2<=0.1)=NaN;
L2 = uint16(L2/L2_max*2^16);
vtkwrite('../../../results/VTKfiles/lambda2N60_t2000.vtk','structured_points','t2000',L2);
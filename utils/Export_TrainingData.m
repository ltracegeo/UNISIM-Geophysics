load('.\Data\TrainingData.mat')


well = read_las_file('.\Data\well_las\INJ003.las');

well.curve_info{3,1} = 'SW_2013'; well.curve_info{3,2} = '_'; well.curve_info{3,3} = 'SW_2013';
well.curve_info{4,1} = 'SW_2024'; well.curve_info{4,2} = '_'; well.curve_info{4,3} = 'SW_2024';
well.curve_info{5,1} = 'Ip2013'; well.curve_info{5,2} = 'kg/(m2.s)'; well.curve_info{5,3} = 'Ip2013';
well.curve_info{6,1} = 'VpVs2013'; well.curve_info{6,2} = '_'; well.curve_info{6,3} = 'VpVs2013';
well.curve_info{7,1} = 'Ip2024'; well.curve_info{7,2} = 'kg/(m2.s)'; well.curve_info{7,3} = 'Ip2024';
well.curve_info{8,1} = 'VpVs2024'; well.curve_info{8,2} = '_'; well.curve_info{8,3} = 'VpVs2024';
well.curve_info{9,1} = ''; well.curve_info{9,2} = ''; well.curve_info{9,3} = '';


depth = well.curves(:,1);
%depth = linspace(depth(1),depth(end),size(TrainingData,1));
depth = 1:150000;

well.step = depth(2)-depth(1);
well.curves = [];
well.curves(:,1) = depth;
well.curves(:,2) = TrainingData(:,5); % por
well.curves(:,3) = TrainingData(:,6); % sw1
well.curves(:,4) = TrainingData(:,7); % sw2
well.curves(:,5) = TrainingData(:,1); % ip1
well.curves(:,6) = TrainingData(:,2); % vp/vs1
well.curves(:,7) = TrainingData(:,3); % ip2
well.curves(:,8) = TrainingData(:,4); % vp/vs2

write_las_file(well,'.\Export\Wells\TraininData_well.las'); 


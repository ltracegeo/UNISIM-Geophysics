function [] = write_segy(dt, seismic, X, Y, t0, name)
%% Inputs: Increments in the x, y and time directions; property data; starting time; file name

Y = Y(1,:,:);
X = X(1,:,:);
tracesx = Y(:)';
tracesy = X(:)';

traces = 1:1:size(seismic,2)*size(seismic,3);
[inlines,xlines] = meshgrid([1:size(seismic,2)],[1:size(seismic,3)]);

seismic = reshape(seismic,size(seismic,1),size(seismic,2)*size(seismic,3));

% Transforms the data into a survey structure
survey = s_convert(seismic,t0,dt);

%% Write Headers
survey=ds_add_header(survey,xlines(:)',{'ds_seqno','n/a','Trace sequence number within line'});
survey=ds_add_header(survey,inlines(:)',{'ffid','n/a','Original Field record number'});
survey=ds_add_header(survey,xlines(:)',{'CDP','n/a','CDP number'});
survey=ds_add_header(survey,tracesy,{'sou_x','m','X coordinate of source'});
survey=ds_add_header(survey,tracesx,{'sou_y','m','Y coordinate of source'});
survey=ds_add_header(survey,tracesy,{'cdp_x','n/a','X-coordinate of CDP'});
survey=ds_add_header(survey,tracesx,{'cdp_y','n/a','Y-coordinate of CDP'});
survey=ds_add_header(survey,inlines(:)',{'iline_no','n/a','In-line number'});
survey=ds_add_header(survey,xlines(:)',{'xline_no','n/a','Cross-line number'});
survey=ds_add_header(survey,1,{'trc_type','n/a','Trace type (1=live,2=dead,3=dummy,...)'});
survey=ds_add_header(survey,1,{'lag','ms','Lag time between shot and recording start'});

write_segy_file(survey, name);

end


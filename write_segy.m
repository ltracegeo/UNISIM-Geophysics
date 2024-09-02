function [] = write_segy(dt, seismic, X, Y, t0, name)
%% Inputs: Increments in the x, y and time directions; property data; starting time; file name

Y = Y(1,:,:);
X = X(1,:,:);
tracesx = Y(:)';
tracesy = X(:)';

traces = 1:1:size(seismic,2)*size(seismic,3);

seismic = reshape(seismic,size(seismic,1),size(seismic,2)*size(seismic,3));

% Transforms the data into a survey structure
survey = s_convert(seismic,t0,dt);

%% Write Headers
survey=ds_add_header(survey,traces,{'ds_seqno','n/a','Trace sequence number within line'});
survey=ds_add_header(survey,tracesy,{'ffid','n/a','Original Field record number'});
survey=ds_add_header(survey,tracesx,{'o_trace_no','n/a','Trace sequence number within original field record'});
survey=ds_add_header(survey,tracesx,{'source','n/a','Energy source point number'});
survey=ds_add_header(survey,tracesx,{'CDP','n/a','CDP number'});
survey=ds_add_header(survey,tracesx,{'seq_cdp','n/a','Trace sequence number within CDP ensemble'});
survey=ds_add_header(survey,tracesy,{'sou_x','m','X coordinate of source'});
survey=ds_add_header(survey,tracesx,{'sou_y','m','Y coordinate of source'});
survey=ds_add_header(survey,tracesy,{'cdp_x','m','X-coordinate of CDP'});
survey=ds_add_header(survey,tracesx,{'cdp_y','m','Y-coordinate of CDP'});
survey=ds_add_header(survey,1,{'trc_type','n/a','Trace type (1=live,2=dead,3=dummy,...)'});

write_segy_file(survey, name);

end


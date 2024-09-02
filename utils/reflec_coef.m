function Rpp=reflec_coef(type, rho_trace,vp_trace, vs_trace, angle)
% function [Rpp,Rzero,G]=reflec_coef(type, rho_trace,vp_trace, vs_trace, angle)
%% COMPUTE THE REFLECTIVITY SERIES TRACE-BY-TRACE FOR NORMAL INCIDENCE AND USING
% SHUEY'S APPROXIMATION AT EACH GATHER LOCATION
% angle-dependent - Shuey's approximation
%
%
% FUNCTION[RC]=REFLEC_COEF(trace1,trace2, trace3, angle)
% type - 1 = NORMAL INCIDENCE;
%        2 = SHUEY LINEAR APPROX.;
%        3 = FULL SHUEY APPROX;
%        4 = AKI-RICHARDS APPROX.
% rho_trace - DENSITY AT k INDEX
% vp_trace - PWAVE-VEL AT k INDEX
% vs_trace - SWAVE-VEL AT k INDEX
%
% angle - 1-DIMENSIONAL ANGLE ARRAY
% rc - OUTPUT REFLECTIVITY SERIES GATHER
%
% MODIFICATIONS
% 12 DEZ 2023: AKI & RICHARDS PERFORMANCE UPDATE
% 03 AUG 2012: diff ADDED TO COMPUTE DIFFERENCES; for CYCLES REMOVED
% 02 AUG 2012: PERFORMANCE UPDATE
% 13 APR 2012: ADDED SHUEY NOT LINEAR; ADDED AKI-RICHARDS;
%              CONSISTENCY CHECK; MODIFICATIONS USING TAPAN'S M-FILE
% 24 MAR 2012: MINOR CHANGES FOR CONSISTENCY. ANGLES IN DEGREES
% 06 MAR 2012: CONSISTENCY CHECK; COMENTARIES ADDED; SIND IS NOW CALCULATED
%              Rzero AND G AS OUTPUTS; LAST LAYER OF RC, Rzero AND G IS NOW
%              SET TO ZERO.
% 27 MAR 2012: INPUT VARIABLES CHANGED. type ADDED
%
% 
% LA - 12 - FEB 2012
%% COMPUTE REFLECTION CONFICIENTS FOR NORMAL INCIDENCE
if type==1
    size_ai=size(rho_trace);
    Rpp=zeros(size_ai(1),1);

    for i=1:size_ai(1)-1
        Rpp(i)=(rho_trace(i+1)-rho_trace(i))/(rho_trace(i+1)+rho_trace(i));
    end
%% COMPUTE REFLECTION COEFFICIENTS FOR SHUEY LINEAR APPROX.
elseif type==2
 
    % ALLOCATE ARRAYS
    den=zeros(length(rho_trace),1);
    Vp=zeros(length(vp_trace),1);
    Vs=zeros(length(vs_trace),1);
    Rzero=zeros(length(vp_trace),1);
    G=zeros(length(vp_trace),1);
    
    Rpp=zeros(length(vp_trace),length(angle));

    % COMPUTE VARIABLES FROM PWAVE/SWAVE/DENSITY MODELS
    Dden=diff(rho_trace);
    DVp=diff(vp_trace);
    DVs=diff(vs_trace);
    Dden(size(Dden,1)+1,1)=0;
    DVp(size(DVp,1)+1,1)=0;
    DVs(size(DVs,1)+1,1)=0;
    for j=1:length(vp_trace)-1
        den(j,1)=(rho_trace(j+1)+rho_trace(j))/2;
        Vp(j,1)=(vp_trace(j+1)+vp_trace(j))/2;
        Vs(j,1)=(vs_trace(j+1)+vs_trace(j))/2;
    end

    Rzero=0.5.*((DVp./Vp)+(Dden./den));
    G=(-2.*Vs.^2.*Dden./(Vp.^2.*den))+(0.5.*(DVp./Vp))-(4*Vs.*DVs./(Vp.^2));
    Rpp=repmat(Rzero,1,size(angle,2))+(G*sind(angle).^2);

    clear Dden DVp DVs den Vp Vs;
    
    Rpp(size(Rpp,1),:)=0;
    Rzero(size(Rpp,1),:)=0;
    G(size(Rpp,1),:)=0;
%%  3 - FULL SHUEY APPROX;
elseif type==3
    
    % ALLOCATE ARRAYS
    den=zeros(length(rho_trace),1);
    Vp=zeros(length(vp_trace),1);
    Vs=zeros(length(vs_trace),1);
    Rzero=zeros(length(vp_trace),1);
    G=zeros(length(vp_trace),1);
    F=zeros(length(vp_trace),1);
    
    Rpp=zeros(length(vp_trace),length(angle));

    % COMPUTE VARIABLES FROM PWAVE/SWAVE/DENSITY MODELS
    Dden=diff(rho_trace);
    DVp=diff(vp_trace);
    DVs=diff(vs_trace);
    Dden(size(Dden,1)+1,1)=0;
    DVp(size(DVp,1)+1,1)=0;
    DVs(size(DVs,1)+1,1)=0;

    for j=1:length(vp_trace)-1
        den(j,1)=(rho_trace(j+1)+rho_trace(j))/2;
        Vp(j,1)=(vp_trace(j+1)+vp_trace(j))/2;
        Vs(j,1)=(vs_trace(j+1)+vs_trace(j))/2;
    end

    % COMPUTE R(0) AND G
    Rzero=0.5.*((DVp./Vp)+(Dden./den));
    G=(-2.*Vs.^2.*Dden./(Vp.^2.*den))+(0.5.*(DVp./Vp))-...
        (4*Vs.*DVs./(Vp.^2));
    F=0.5.*(DVp./Vp);
    Rpp=repmat(Rzero,1,size(angle,2))+(G*sind(angle).^2)+(F*(tand(angle).^2-sind(angle).^2));

    Rpp(size(Rpp,1),:)=0;
    clear Dden;
    clear DVp;
    clear DVs;
    clear den;
    clear Vp;
    clear Vs;
    
%% 4 - AKI-RICHARDS APPROX.    
elseif type==4
        
    % ALLOCATE ARRAYS
    den=zeros(length(rho_trace),1);
    Vp=zeros(length(vp_trace),1);
    Vs=zeros(length(vs_trace),1);
    ct=cosd(angle);
    
    % COMPUTE VARIABLES FROM PWAVE/SWAVE/DENSITY MODELS
    Dden=diff(rho_trace);
    DVp=diff(vp_trace);
    DVs=diff(vs_trace);
    Dden(size(Dden,1)+1,1)=0;
    DVp(size(DVp,1)+1,1)=0;
    DVs(size(DVs,1)+1,1)=0;
   
    for j=1:length(vp_trace)-1
        den(j,1)=(rho_trace(j+1)+rho_trace(j))/2;
        Vp(j,1)=(vp_trace(j+1)+vp_trace(j))/2;
        Vs(j,1)=(vs_trace(j+1)+vs_trace(j))/2;
    end
    p=repmat(sind(angle),size(vp_trace,1),1)./repmat(vp_trace,1,size(angle,2));  
     Rpp=(0.5.*(1-(4.*p.^2.*repmat(Vs,1,size(angle,2)).^2)).*repmat(Dden,1,size(angle,2))./repmat(den,1,size(angle,2)))+...
        (repmat(DVp,1,size(angle,2))./(2.*repmat(ct,size(vp_trace,1),1).^2.*repmat(Vp,1,size(angle,2))))-...
       (4.*p.^2.*repmat(Vs,1,size(angle,2)).*repmat(DVs,1,size(angle,2)));
     %Rpp=(0.5.*(1-(4.*p.^2.*repmat(vs_trace,1,size(angle,2)).^2)).*repmat(Dden,1,size(angle,2))./repmat(den,1,size(angle,2)))+...
     %   (repmat(DVp,1,size(angle,2))./(2.*repmat(ct,size(vp_trace,1),1).^2.*repmat(Vp,1,size(angle,2))))-...
     %  (4.*p.^2.*repmat(vs_trace,1,size(angle,2)).*repmat(DVs,1,size(angle,2)));
   
    Rpp(size(Rpp,1),:)=0;
    Rzero=0;
    G=0;
    clear Dden;
    clear DVp;
    clear DVs;
    %% 5 - FATTI APPROX.    
elseif type==5
        
    % ALLOCATE ARRAYS
    den=zeros(length(rho_trace),1);
    Vp=zeros(length(vp_trace),1);
    Vs=zeros(length(vs_trace),1);
    ct=cosd(angle);
    
    % COMPUTE VARIABLES FROM PWAVE/SWAVE/DENSITY MODELS
    Dden=diff(rho_trace);
    DVp=diff(vp_trace);
    DVs=diff(vs_trace);
    Dden(size(Dden,1)+1,1)=0;
    DVp(size(DVp,1)+1,1)=0;
    DVs(size(DVs,1)+1,1)=0;
   
    for j=1:length(vp_trace)-1
        den(j,1)=(rho_trace(j+1)+rho_trace(j))/2;
        Vp(j,1)=(vp_trace(j+1)+vp_trace(j))/2;
        Vs(j,1)=(vs_trace(j+1)+vs_trace(j))/2;
    end
    p=repmat(sind(angle),size(vp_trace,1),1)./repmat(vp_trace,1,size(angle,2));  
    Rpp=(0.5.*(tand(angle).^2-(4.*p.^2.*repmat(Vs,1,size(angle,2)).^2)).*repmat(Dden,1,size(angle,2))./repmat(den,1,size(angle,2)))+...
    (repmat(DVp,1,size(angle,2))./(2.*repmat(ct,size(vp_trace,1),1).^2.*repmat(Vp,1,size(angle,2))))-...
(4.*p.^2.*repmat(Vs,1,size(angle,2)).*repmat(DVs,1,size(angle,2)));
 

    Rpp(size(Rpp,1),:)=0;
    Rzero=0;
    G=0;
    clear Dden;
    clear DVp;
    clear DVs;
%% 6 - FOSTER APPROX.    
elseif type==6
        
    % ALLOCATE ARRAYS    
    Vp=zeros(length(vp_trace),1);
    Vs=zeros(length(vs_trace),1);
    ct=cosd(angle);
    st=sind(angle);
    
    % COMPUTE VARIABLES FROM PWAVE/SWAVE/DENSITY MODELS
    DVp=diff(vp_trace);
    DVs=diff(vs_trace);
    DVp(size(DVp,1)+1,1)=0;
    DVs(size(DVs,1)+1,1)=0;
   
    for j=1:length(vp_trace)-1
        Vp(j,1)=(vp_trace(j+1)+vp_trace(j))/2;
        Vs(j,1)=(vs_trace(j+1)+vs_trace(j))/2;
    end
    
    Rpp = (0.5*repmat(DVp,1,size(angle,2))./repmat(Vp,1,size(angle,2))).*repmat(ct,size(vp_trace,1),1).^2  + ...        
            (repmat(DVs,1,size(angle,2))./repmat(Vs,1,size(angle,2))).*repmat(st,size(vp_trace,1),1).^2;



        
   
    
    Rpp(size(Rpp,1),:)=0;
    Rzero=0;
    G=0;
    clear Dden;
    clear DVp;
    clear DVs;
end
end

    
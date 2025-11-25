% function Matrix = Heat_MultipleScans(location,state,style,scan_pattern)
% 
% SS = 600; %scan speed 400 - 1500
% LP = 100; %Laser Power 100 - 300
% eeta = 0.3; %Absorptivity around 0.4
% r_b = 0.06; %Laser Radius
% H = 0.1; % Hatch distance < 2 * r_b
% 
% scans = 10;
% 
% if(isnan(state.time))
%     Matrix = nan(1, length(location.x));
% 
% else
%     flag = 1;
%     time_ind = ceil(state.time/0.005); %scans (time/scan time for each line)
%     if(time_ind == 0)
%         time_ind = 1;
%     end
% 
%     if(time_ind > scans)
%         time_ind = scans;
%         flag = 0;
%     end
% 
%     x_centre =  0.1 + (H*(scan_pattern(time_ind)));
% 
% 
%     x_diff2 = (location.x - x_centre).^2;
% 
% 
%     r_b2 = (r_b).^2;
% 
% 
%     Matrix = 6*sqrt(3)*LP*eeta*exp(-3*(x_diff2)./r_b2)./(pi*r_b);%Goldak's Heat source
% 
%     Matrix = Matrix.*flag; 
% 
% 
% end
% end

function Matrix = Heat_MultipleScans(location,state,params,style,scan_pattern)

    SS=params.SS; % scan speed
    LP=params.LP;   % laser power
    eeta=params.eeta;   % absorptivity
    r_b=params.r_b;  % beam radius
    H=params.H;   % hatch distance

    scans = length(scan_pattern);

    if isnan(state.time)
        Matrix = nan(1, numel(location.x));
        return
    end

    if style == "simultaneous"
        % disp("Printing simultaneously")

        flag = 1;
        time_ind = ceil(state.time/0.005);
        if time_ind == 0, time_ind = 1; end
        if time_ind > scans
            time_ind = scans;
            flag = 0;
        end
        x_centres = scan_pattern;
        x = location.x(:);
        rb2 = r_b^2;
    
        xdiff2 = (x - x_centres).^2; 
    
        G = exp(-3 * xdiff2 ./ rb2);
        amp = 6*sqrt(3)*LP*eeta / (pi*r_b);
        Matrix_col = amp * sum(G, 2);
    
        Matrix = (Matrix_col.' ) * flag;

    elseif style == "sequential"

        % disp("Printing sequentially")
        flag = 1;
        time_ind = ceil(state.time/0.005); %scans (time/scan time for each line)
        if(time_ind == 0)
            time_ind = 1;
        end
    
        if(time_ind > scans)
            time_ind = scans;
            flag = 0;
        end
    
        x_centre =  (scan_pattern(time_ind));
    
        x_diff2 = (location.x - x_centre).^2;
    
        r_b2 = (r_b).^2;    
        Matrix = 6*sqrt(3)*LP*eeta*exp(-3*(x_diff2)./r_b2)./(pi*r_b);%Goldak's Heat source
        Matrix = Matrix.*flag;
       
    else
        disp('Unknown styles')
    end

end
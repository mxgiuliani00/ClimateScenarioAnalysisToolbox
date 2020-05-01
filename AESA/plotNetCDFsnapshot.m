function Xt = plotNetCDFsnapshot(x, dstart, dend, td, var)

% Xt = plotNetCDFsnapshot(x, dstart, dend, td, var)
%
% This function produces a figure showing a map of the variable stored in
% the matlab structure x (created via extractNetCDFdata function) for a
% given day (td) in the time horizon (dstart-dend), which is associated to
% the scenario of x.
%
% Input:    - x = matlab structure obtained from a NetCDF file using the
%                   extractNetCDFdata function 
%           - dstart = first day of the time horizon associated to the
%           scenario of x (vector [yyyy, mm, dd])
%           - dend = last day of the time horizon associated to the
%           scenario of x (vector [yyyy, mm, dd])
%           - td = specific date for taking the snapshot (vector [yyyy, mm, dd])
%           - var = type of variable for setting colormap (string).
%                   the function currently works only for precipitation (pr) 
%                   and temperature (tas or ts)
% Output:   - Xt = visualized 2D snapshot matrix 
%
% Last Update: MatteoG, 11/12/2015

% create date vectors and remove Feb. 29
dn_hist = datenum(dstart(1),dstart(2),dstart(3)):datenum(dend(1),dend(2),dend(3));
dv_hist = datevec(dn_hist);
id29feb = dv_hist(:,2).*dv_hist(:,3) == 2*29;
dv_hist_nl = dv_hist(~id29feb,1:3) ;
dn_hist_nl = datenum(dv_hist_nl(:,1),dv_hist_nl(:,2),dv_hist_nl(:,3));

% identify index of selected day (td)
tdn = datenum(td(1),td(2),td(3)) ;
if tdn < dn_hist_nl(1) || tdn > dn_hist_nl(end)
    error('the selected date is outside the time horizon of this scenario')
else
    idx_td = find(dn_hist_nl == tdn) ;
    Xt = x.value(:,:,idx_td);
    % create figure
    if strcmp(var,'tas') || strcmp(var,'ts')
        colormap(jet);
        pcolor( x.lon, x.lat, Xt );
        shading interp; colorbar;
    elseif strcmp(var,'pr')
        colormap(flipud(winter));
        pcolor( x.lon, x.lat, Xt );
        shading interp; colorbar;
    else
        error('the selected variable is neither precipitation (pr) or temperature (tas,ts)')
    end
end

end

% Copyright 2015 Matteo Giuliani, Politecnico di Milano
% M. Giuliani: matteo.giuliani@polimi.it - http://giuliani.faculty.polimi.it
% Remove fit information from DMS glass variables, so that they will be fit
% again.
pareval('Glass','if(isfield((),''muX'')) () = rmfield((),fields); end')
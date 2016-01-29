function Label = PlotLabel(FieldName)
% Mapping from FieldName to label strings for commonly plotted fields

switch(FieldName)
    case 'CorrX'
        Label = 'Magnetization (emu/cc)';
    case 'RawCorrX'
        Label = 'Moment (uemu)';
    case 'Field'
        Label = 'Applied Field (Oe)';
    case 'Kerr'
        Label = 'Normalized Kerr Signal';
    case 'Kerr_x'
        Label = 'Kerr Signal (Blue)';
    case 'Kerr_y'
        Label = 'Kerr Signal (Red)';
    otherwise
        Label = FieldName;
end

end


function f_print_report(fid, titleText, R)
% ==================================================================
% f_print_report(fid, titleText, R)
%
% Writes a report cell array R = {label, value, note; ...} (Report_P1 ... P5)
% to the file identifier fid (use fid = 1 for the command window).
% ==================================================================
fprintf(fid, '\n%s\n', titleText);
for i = 1:size(R,1)
    v = R{i,2}; note = R{i,3};
    if isnumeric(v) || islogical(v)
        vs = mat2str(double(v), 6);
    elseif ischar(v)
        vs = v;
    else
        vs = '';
    end
    if ~ischar(note), note = ''; end
    fprintf(fid, '  %-48s %-42s %s\n', R{i,1}, vs, note);
end
end

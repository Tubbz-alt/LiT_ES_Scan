function disp_prog(current,total,del)

percent = 100*current/total;
prog_string = [num2str(percent,'%05.2f') '%%'];
if del; fprintf('\b\b\b\b\b\b'); else fprintf('Progress: '); end;
fprintf(prog_string);
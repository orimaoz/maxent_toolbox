% prints a pass/fail status (for the unit testing)
function printStatus(str_status)

switch lower(str_status)
    case 'pass'
        cprintf([0,0.5,0],[str_status '\n']);       
    case 'partial'
        cprintf([0.6 0.6 0],[upper(str_status) '\n']);                       
    case 'fail'
        cprintf('err',[upper(str_status) '\n']);               
    otherwise
        cprintf([str_status '\n']);
end

end

function progressbar
    for i = 0:100
        progressbari(i,100);
        pause(.1);
    end

function progressbari(CurrentVal,MaxVal)
if CurrentVal == 0;
    fprintf(1,'[          ]');
else
    progress = floor(CurrentVal*10 / MaxVal); % 0-10
    
    fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b'); % Erase bar
    for j = 1:progress
        fprintf(1,'-'); % Print 1-10 '-'
    end
    for k = 1:(10-progress)
        fprintf(1,' ');
    end
    fprintf(1,']');
    
    if CurrentVal == MaxVal
       fprintf('\n') 
    end
end
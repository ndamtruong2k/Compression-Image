function [ times ] = time( ending )
%{
By Jacob Thomas
Gives time taken to complete a given operation must be started using tic 
and ending must use toc.
 ex
start=tic;
---- Code being timed ---
ending=toc(start);
[times]=time(ending)
-----------------------------
Outputs
times=
-- hr -- min --.-- sec
%}

times=[num2str(fix(ending/3600)),' hr '...
    ,num2str(fix((ending-fix(ending/3600)*3600)/60)),' min '...
    ,num2str(ending-fix((ending-fix(ending/3600)*3600)/60)*60-...
    fix(ending/3600)*3600),' sec'];
end


function yref = myRef_tracking(t, ts)
% helper function to define tracking references; this function should be
% customized to the desired demonstration, or an equivalent function may be
% designed.
Nstep = 1*60/ts;

if ceil(t/ts) < Nstep
    yref = 40.0;
elseif ceil(t/ts) < 2*Nstep
    yref = 38.0;
elseif ceil(t/ts) < 3*Nstep
    yref = 42.0;
elseif ceil(t/ts) < 4*Nstep
    yref = 36.0;
else
    yref = 40.0;
end

end
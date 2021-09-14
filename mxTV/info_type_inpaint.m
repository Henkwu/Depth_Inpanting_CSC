function info = info_type_inpaint(type,N,k,eps_rel_k,time)
%
% Type-like definitions for TVinpaint to be used with mxTV software
%
    
    info_type = struct('STATUS',{'EPSILON-OPTIMAL-SOLUTION',...
        'NOT-EPSILON-OPTIMAL-SOLUTION',...
        'MAXIMUM-NUMBER-OF-ITERATIONS-EXCEEDED'},...
        'NINPAINT',{N},'ITERATIONS_K',{k},...
        'EPS_REL_K',{eps_rel_k},...
        'TIME',{time});

    info = info_type(type);
function any_error = install_windows(fftwpath,libfile)
%INSTALL_WINDOWS Script used to install mxTV on the Windows platform.
%
% Compiles and links the mex files for the mxTV package.
%
% Usage: install_windows(<fftwpath>,<libfile>), where
% <fftwpath> is a string with the path to where fftw was placed, and 
% <libfile> is the path to an fftw import library file (*.lib)
% that can be used with the selected compiler.
%
% The file externlib/Lcc/libfftw3-3.lib is an import library file
% for the Lcc compiler
%
% The file externlib/msdk/libfftw3-3.lib is an import library file
% for the Microsoft Software Development Kit (SDK)
%
% Tested with the Lcc compiler bundled with older versions of Matlab.
% and Microsoft Software Development Kit (SDK)
%
% See readme.txt for further instructions.
%
% J. Dahl^1, P.C. Hansen^2, S.H. Jensen^1 & T.L. Jensen^1
% CSI project: (1) Aalborg University, (2) Technical University of Denmark
% April 28, 2009.
%
% If you want to be able to break the execution of the programs, try to
% set CTRLCBREAK = 1, which uses a non-documented MATLAB API.
% If you do not have libut.lib in your Matlab dir, try CTRLCBREAK = 2.
% Default CTRLCBREAK=0.
%

CTRLCBREAK=0;

if CTRLCBREAK==0
    sbreak = '';
elseif CTRLCBREAK==1
    sbreak = ['-DLIBUT -L' matlabroot '\extern\lib\win32\lcc -llibut'];
elseif CTRLCBREAK==2
    sbreak = ['-DLIBUT -Lexternlib -llibut'];
else
    error('Not a valid option for CTRLCBREAK')
end


ext = mexext;
if strcmp(ext(end-1:end),'64')
	arg = '-largeArrayDims';
else
	arg = '';
end

any_error = false;

% Denosing.
try
    cs = sprintf('mex %s %s c/tv_denoise.c c/tools.c c/tv_core.c',sbreak,arg);
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVdenoise because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end

% Inpainting.
try
    cs = sprintf('mex %s %s c/tv_inpaint.c c/tools.c c/tv_core.c',sbreak,arg);
    eval(cs)
catch
    any_error= true;    
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVinpaint because the compilation failed.')
    disp('Follow the above instructions to locate the problem.')
end


% Deblurring.
deblur_error= false;
try
    cs = sprintf('mex -I%s %s c/mxtrp.c c/tools.c %s',fftwpath,arg,libfile);
    eval(cs)	

    cs = sprintf('mex -DFFTW3 -I%s %s %s c/tv_deblur.c c/tv_core.c c/tools.c %s',fftwpath,sbreak,arg,libfile);
    eval(cs)
  
    cs = sprintf('mex -DFFTW3 -I%s %s %s c/tv_deblur_rr.c c/tv_core.c c/tools.c %s',fftwpath,sbreak,arg,libfile);
    eval(cs)
catch
    any_error= true;
    disp('-------------------------------------------------------------------')
    disp('You will not be able to use TVdeblur because the compilation failed,')
    disp('or linking failed. To locate the problem, follow the above') 
    disp('instructions or see the readme.txt file.')
    disp('Ignore this error if you do not need TVdeblur.')
    deblur_error = true;
end

%We will just end by trying out if fftw has the correct setup
if deblur_error == false
	try	
        xtemp=TVdeblur(1,1,0.1);
        clear xtemp
    catch
        any_error = true;
        disp('----------------------------------------------------------------')
        disp('You will not be able to use TVdeblur because the setup of fftw3')
        disp('is wrong. Follow the instructions in readme.txt')
        disp('to solve the problem if you would like to use TVdeblur.')
        disp('Ignore this error if you do not need TVdeblur.')
	end
end

if any_error == false && deblur_error == false
    disp('Install completed successfully.')
else
    disp('Installation did NOT complete successfully.')
end

%
%  Simple Optimizer that connects to CgMx
%
% Usage:
%    optimizer -caseName=[cyl|block|blockWidth|lens] -plotOption=[0|1] -tf=<f> -probeType=[point|transmission] ...
%              -method=[fake,fminsearch] -gridFactor=<i> -plotGrid=[0|1] -plotSolution=[0|1] ...
%              -objective=[minimizeReflection|targetTransmission] -targetFile=<>
% 
% caseName:
%     block        : scattering from a dielectric block, change epsilon
%     blockWidth   : scattering from a dielectric block, change the width 
% -plotGrid = 1 : plot the grid afet each optimizer step
% -plotSolution = 1 : plot the solution after each optimizer step
% 
% Examples
%
% 
function optimizer(varargin)

 fontSize=14;  lineWidth=2;  markerSize=5;   % for plots
 
 % Define some global variables to avoid pass so many args to runMaxwell
 globalDeclarations

 method = 'fake'; 

 caseName = 'cyl';   % case to run 
 plotOption=1;
 plotGrid=1;
 plotSolution=1;
 infoLevel=1;
  
 objective='minimizeReflection'; 
 targetFile='none'; % data file for target 

 pointProbe=0; transmissionProbe=1; % probe types 
 probeType='point';
 tFinal=1;         % final time 
 gridFactor=4;     % defines grid resolution 1,2,4,8
 kx=1;             % wave number 
 blockWidth=.5;    % default width of dielectric block
 eps1=4.;          % default block epsilon
 tolFun=.1; tolX=.1; % tolerences for fminserach
 x0Default=-1234567.; 
 x0 = x0Default;     % initial guess  set to a real number to change the default initial guess

  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};
    caseName      = getString( line,'-caseName',caseName );
    probeType     = getString( line,'-probeType',probeType );
    method        = getString( line,'-method',method );
    objective     = getString( line,'-objective',objective );
    targetFile    = getString( line,'-targetFile',targetFile );
    tFinal        =   getReal( line,'-tf',tFinal );
    kx            =   getReal( line,'-kx',kx );
    x0            =   getReal( line,'-x0',x0 );
    eps1          =   getReal( line,'-eps1',eps1 );
    tolFun        =   getReal( line,'-tolFun',tolFun );
    tolX          =   getReal( line,'-tolX',tolX );
    blockWidth    =   getReal( line,'-blockWidth',blockWidth );
    infoLevel     =    getInt( line,'-infoLevel',infoLevel );
    plotOption    =    getInt( line,'-plotOption',plotOption );
    plotGrid      =    getInt( line,'-plotGrid',plotGrid );
    plotSolution  =    getInt( line,'-plotSolution',plotSolution );
    gridFactor    =    getInt( line,'-gridFactor',gridFactor );
  end

  fprintf('>>> optimizer: caseName=%s, method=%s, \n objective=%s (targetFile=%s, tolFun=%g, tolX=%g), infoLevel=%d, \n plotOption=%d, plotGrid=%d, plotSolution=%d, gridFactor=%d, tFinal=%9.3e, probeType=%s, \n',...
                caseName,method,objective,targetFile,tolFun,tolX,infoLevel,plotOption,plotGrid,plotSolution,gridFactor,tFinal,probeType);
  fprintf('             : kx=%g, blockWidth=%g, eps1=%g\n',kx,blockWidth,eps1);


  iteration=0;      % keeps track of how many times runMaxwell is called. 
  maxIterations=10; 

  % ----------------------------------------------------
  %   Define an objective function for minimization
  % ----------------------------------------------------
  function ff = cgmxFunction( xx )
    % fprintf(' objectiveFunction: x=%g\n',xx(1)); 
    % fprintf(' objectiveFunction: size(xx,2)=%d\n',size(xx,2));

    if( strcmp(objective,'targetTransmission') )

      % -- LENS: for now just adjust the right control point: 
      par(1)=kx;
      par(2)=eps1; 

      if( size(xx,2)==1 )
        % shift right control point
        dxLeft =-.2; dxRight=.2 + xx(1); 
      else
        % shift left and right control points 
        dxLeft =-.2 + xx(1); dxRight=.2 + xx(2); 
      end

      par(3)=dxLeft;    % shift left control point
      par(4)=dxRight;   % shift right control point
	
    elseif( strcmp(objective,'minimizeReflection') )
      
      if( strcmp(caseName,'blockWidth') )
        blockWidth=xx(1); 
      else
        eps1=xx(1);
      end 
      par(1)=kx;
      par(2)=eps1;
      par(3)=blockWidth; 

    else

      fprintf('optimizer: ERROR: unknown objective =[%s]\n',objective);
      pause; pause;
    end;


      
    [ values ]  = runMaxwell( caseName,tFinal,probeType,gridFactor,infoLevel,plotOption, par );
    ff=values(1); % reflection coefficient 

  end

  % ff = @(x) x(1)^2;
  ff = @(x) cgmxFunction( x ); % for some reason we need to call cgMxFunction this way

  if( strncmp(method,'fminsearch',length('fminsearch')) )

    % ---- Find minium using fminsearch: ----

    options = optimset('Display','iter','TolFun',tolFun, 'TolX',tolX,'PlotFcns',@optimplotfval);

    % ------ INITIAL GUESS FOR FMINSEARCH -----
    if( strcmp(objective,'targetTransmission') )

      % --- Adjust the shape of the lens: 

      if( x0 == x0Default  )
        x0 = [ -.05, .05  ]; % initial guess for dxRight (offset from 'exact') 
      else
        x0 = [ x0, -x0 ];  % user specified initial guess 
      end; 	

    elseif( strcmp(objective,'minimizeReflection') )
      
      if( strcmp(caseName,'blockWidth') )
        x0=[ blockWidth ]; %  blockWidth
      else
        x0=[ 3.]; %  3. ,  8. initial guess
      end;
    else

      fprintf('optimizer: ERROR: unknown objective =[%s]\n',objective);
      pause; pause;
    end;
    
    % fval = myObjectiveFunction( x0 );
    % fprintf('myObjectiveFunction: x0=%g, f=%g\n',x0(1),fval);
    
    fprintf('call fminsearch...\n');
    figure(2); 

    % [x,fval] = fminsearch(myObjectiveFunction,x0,options);
    [x,fval] = fminsearch(ff,x0,options);

    if( size(x,2)==1 )
      fprintf('...DONE fminsearch: x=%g, fval=%g\n',x(1),fval);
    else
      fprintf('...DONE fminsearch: x=[%g,%g], fval=%g\n',x(1),x(2),fval);
    end; 
    if( plotOption > 0  )
      grid on; set(gca,'FontSize',fontSize);
      plotFileName=sprintf('%sFminsearchConvergence.eps',caseName);
      fprintf('optimizer: save plot file=[%s]\n',plotFileName);
      print('-depsc2',plotFileName); % save as an eps file
    end

  else

    % -------------------- FAKE OPTIMIZATION LOOP  ------------------------
    for iter=1:maxIterations
    
      if( strcmp(caseName,'blockWidth') ) 
        kx=1;
        eps1=4;
	blockWidthNew = blockWidth+.1*iter;
        par(1)=kx;
        par(2)=eps1; 
        par(3)=blockWidthNew

      elseif( strcmp(caseName,'block') ) 
        kx=1;
        eps1= 2+ .5*iter;
        par(1)=kx;
        par(2)=eps1; 
        par(3)=blockWidth

      elseif( strcmp(caseName,'lens') ) 

        % NOTE: for now we just shift one control point at the center of the left and right sides.

        par(1)=kx;
        par(2)=eps1; 
        dxLeft = -.2 + .025*(iter-1); dxRight=.2 - .025*(iter-1); 
        par(3)=dxLeft;    % shift left control point
	par(4)=dxRight;   % shift right control point 

      else
  
        eps1=4; 
        kx = 2 + iter; % incident wave number 
        par(1)=kx;
        par(2)=eps1; 
        par(3)=blockWidth
      end;
      
      [ values ]  = runMaxwell( caseName,tFinal,probeType,gridFactor,infoLevel,plotOption, par );
  
      fprintf('optimizer: iter=%d: return-values=[%g,%g]\n',iter,values(1),values(2));

      if( 1==1 || ( plotGrid==0 && plotSolution==0) )
         pause
      end    
    end; % end for iter
  end; 

fprintf('done\n'); 

end

% --- Utility functions ---





% Function getReal: read a command line argument for a real variable
function [ val ] = getReal( line,name,val)
 % fprintf('getReal: val=%g line=[%s] name=[%s]\n',val,line,name);
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%e',name)); 
   % fprintf('getReal: scan for val=%g\n',val);
 end
end

% Function getInt: read a command line argument for an integer variable
function [ val ] = getInt( line,name,val)
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%d',name)); 
 end
end

% Function getString: read a command line argument for a string variable
function [ val ] = getString( line,name,val)
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%s',name)); 
 end
end

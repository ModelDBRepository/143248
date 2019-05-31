% This function computes one of the 4 profiles described in the paper for
% EC spikes
% Luisa Castro, FCUP

function[y]=intensity_func(x,profile,fmin,fmax,B_mod,A_mod,f,miu,sig,fini,T)

Tstart = 200;            %[ms]
Tstop  = T-Tstart;       %[ms]
[ax,bx]=size(x);

if profile=='G'
    % Gaussian profile with f period sinusoidal component  
    y=fmin*(1+A_mod*cos(2*pi*f*x+fini)).*(1+B_mod*exp((-((x-miu)/sig).^2)/2));
elseif profile=='S'
    % Square Wave profile
    for i=1:bx
        if x(i)>Tstart && x(i)<Tstop
            y(i) = fmax;
        else
            y(i) = fmin;
        end
   end
elseif profile=='P'
    % Pyramid Wave profile
    for i=1:bx
        if x(i)>=0 && x(i)<T/2
            y(i)=(fmax-1)*(x(i)-0)/(T/2-0)+1;
        else
            y(i)=-(fmax-1)*(x(i)-T/2)/(T-T/2)+fmax;
        end
    end
elseif profile=='R'
    % Ramping profile
    Tstart = -300;              % estimated value for spikes to start aprox. at 200ms for T=1000ms 
    for i=1:bx
        if x(i)>=Tstart && x(i)<Tstop
            y(i)=(fmax-fmin)*(x(i)-Tstart)/(Tstop-Tstart)+fmin;
        else
            y(i) = fmin;
        end
    end
end

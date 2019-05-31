%Generates realizations of a Non Homogeneous Poisson Process
%(ie, variable lambda parameter) with a certain probability/intensity
%function (described in intensity_func)
% Luisa Castro, FCUP.


function[trainin]=PoissonProc_Generator(T,dt,profile,fmin,fmax,B_mod,A_mod,f,miu,sig,fini)

trainin=[];             % vector of firing times

%Implementation of the algorithm of Ross in book Simulation (page 85)

k=50;                   % number of segments to partition T

vec=[];                 % setting the k intervals to major with lambdau

lambdau=[];             % variable to store the maximum of each interval

S=[];                   % auxiliary variable

% Constructing the partition
for i=1:k
    vec(i,:)=[(i-1)*(T/k)+dt:dt:(i)*(T/k)];
end


% Computing the maximum for each interval (partition)
for i=1:k
    lambdau(i)=max(intensity_func(vec(i,:),profile,fmin,fmax,B_mod,A_mod,f,miu,sig,fini,T));
end

% Computing one realization of the Non Homogeneus Poisson Process
t=0;
J=1;
I=0;
while  t<=T
    u1=rand;
    X=-log(u1)/lambdau(J);
    if t+X>max(vec(J,:));
        if J==k
            break
        else
            X=(X-max(vec(J,:))+t)*lambdau(J)/lambdau(J+1);
            t=max(vec(J,:));
            J=J+1;
        end
    else
        t=t+X;
        u2=rand;
        if u2<=intensity_func(t,profile,fmin,fmax,B_mod, A_mod,f,miu,sig,fini,T)/lambdau(J);
            I=I+1;
            S(I)=t;
        end
    end
end

trainin=S;
trainin=trainin(trainin<T); %Delete times >T

% % Histogram of EC Firing Times 
% figure
% hist(trainin,T)
% xlabel('Time [ms]')


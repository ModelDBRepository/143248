function new_P = unwindphases( T, P, t0, p0 )
%UNWINDPHASES       unwind a set of cyclic phases
%
%   Description:
%   T is a vector of times
%   P is a vector of phases
%   t0 and p0 are used to set the border line between points ribbons:
%       t0 is the point where the border crosses the phase axis
%       p0 is the point where the border crosses the time axis
%
%  2pi ---------------------...
%     |  ..     \    ..
%   p0|\  :.      \    :.
%     |  \  ::.     \    .::
%     |    \  .::.    \     .:
%     |   .  \  ..      \
%     |    :.  \  :..     \
%   0  ----------------------...
%     0         t0
%
%   Author: Paulo Aguiar, pauloaguiar@fc.up.pt

m = - p0/t0;   % slope

for i = 1:length(T)
    
    % set border equation, of type y = m*x + b
    k = fix( T(i)/t0 );   % border number
    b = (k * 2*pi) + p0;  % border line parameter
    offset = 2.0*pi * ( 1 - k );
    if P(i) < m * T(i) + b      % below border
        new_P(i) = P(i) + offset;
    else                        % above border
        new_P(i) = P(i) + offset - 2.0*pi;
    end
    
end

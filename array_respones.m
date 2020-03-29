function [ a ]=array_respones( azimuth, N, d, lamada ) 
a=[ ];
for i=1 : length( azimuth )
% a=[a (exp(1i*[0:N-1]*2*pi*d*sin(azimuth(i))/lamada)).'];
% a=[ a (sqrt( 1/N ) * exp( 1i * [ 0 : N-1 ] * 2 * pi * d * sin( azimuth(i) )/lamada)).' ];
a=[ a (sqrt( 1/N ) * exp( 1i * [ -(N-1)/2 : (N-1)/2 ] * 2 * pi * d * sin( azimuth(i) )/lamada)).' ];
end

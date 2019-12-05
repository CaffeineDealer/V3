function [m,p,r]=polarplot(V,SponRate,LineStyle)


if nargin<3
    LineStyle = 'b-';
end

V = make_row(V);
n=length(V);
x=[0:2*pi/n:2*pi];           % x is: pi/4, pi/2, 3pi/4...

if 4==n
    x=x+pi/4;
    [m,p,r]=rayleighz2(V);   % mean, p value and r is std dev.
else
    [m,p,r]=rayleighz(V);
end;


V(n+1)=V(1);
% x
% V
Handle = polar(x,V,LineStyle);                   % will make the graph (Theta is x) and radius is V
set(Handle, 'Linewidth', 2);

hold on;
a=[0 r];
% a=[0 r*mean(max(V))];
% a=[0 r*max(V)];
b=[0 m];
Handle = polar(b,a,LineStyle);
set(Handle, 'LineWidth', 3);

if nargin>=2 && ~isempty(SponRate)
    thetas = 0:0.1:2*pi;
    circlex = SponRate*cos(thetas);
    circley = SponRate*sin(thetas);
    plot(circlex, circley,'k--','LineWidth', 2); hold on;
end
% legend('First','Second','Third');
hold off;


return;

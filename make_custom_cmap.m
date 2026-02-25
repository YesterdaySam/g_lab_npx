function [customMap] = make_custom_cmap(col1,col2,n)

R = linspace(col1(1),col2(1),n);
G = linspace(col1(2),col2(2),n);
B = linspace(col1(3),col2(3),n);

customMap = [R', G', B'];
end
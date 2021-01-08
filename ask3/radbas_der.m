%derivative of radbas, gives the dF(p), where p all the input points
function [d] = radbas_der(p)
d = -2*p.*exp(-p.*p);

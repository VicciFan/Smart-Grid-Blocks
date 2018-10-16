function V = flat_start(n_nodes,n_phases)

V = repmat(exp(-1i*(0:1:(n_phases-1))'*2*pi/n_phases),n_nodes,1);

end
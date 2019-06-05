function [w] = opt_beamformer (h)

[A,B] = eig(h*h');
[M,I] = max(diag(B));
w = A(:,I);

end
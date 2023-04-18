function coh = normalizeCSD(csd)
% given a cross-spectral density, covariance, or precision, this function
% normalizes the matrix using the diagonals
% ========================================
% INPUT
% csd: cross-spectral density, precision or covariance
% ========================================
% coh: normalized matrix

coh = (abs(csd).^2) ./ (diag(csd) * diag(csd)');
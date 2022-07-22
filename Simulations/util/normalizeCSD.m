function coh = normalizeCSD(csd)

coh = (abs(csd).^2) ./ (diag(csd) * diag(csd)');
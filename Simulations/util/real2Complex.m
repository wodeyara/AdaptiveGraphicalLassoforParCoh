function covMat = real2Complex(covMat, flag)
% Function takes in the covariance matrix in either real or complex form
% and converts from one to the other.
% ====================================
% INPUT: 
% covMat - augmented covariance matrix
% flag - 0 is covariance matrix is real or 1 if complex
% ====================================
% OUTPUT
% covMat - cross-spectral density / covariance matrix
% For more details see Chapter 2.1 of Schreier and Scharf 2010 

    n = length(covMat)/2;
    T = [eye(n), 1i*eye(n); eye(n), -1i*eye(n)];
    
    if flag == 1 %implies the covMat is complex
        covMat = .25*T' * covMat *T;
    else %covMat is real
        covMat = T*covMat *T';
    end
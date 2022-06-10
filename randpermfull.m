function R = randpermfull(N)
% RANDPERMFULL - Random complete or full permutation (derangement)
%   RANDPERMFULL(N) is a random complete or full permutation of the
%   integers from 1 to N; each integer is randomly assigned to an index
%   position that is different from its value.  In other words R =
%   RANDPERMFULL(N) creates a row vector R, such that R(j) ~= j for all
%   j=1:N.  This type of permutation is called a derangement. 
%   N should be larger than 1.
%
%   Examples:
%     randpermfull(7) % might be [2 4 5 6 1 7 3], but not [2 1 3 6 4 7  5]
%                     % since 3 cannot be in index position 3.
%   
%     % Scramble a sentence
%     Words = {'can','you','make','sense','of','this','scrambled','sentence'} ;
%     R = randpermfull(numel(Words)) ;
%     disp([sprintf('%s ',Words{R}) '?']) 
%
%   Notes:
%   - RANDPERMFULL calls RANDPERM and therefore changes the state of the
%     random number generator.
%   - N should be larger than 3 to obtain a real random full permutation;
%     when N is 2, this function always [2 1]
%   - Reference: http://en.wikipedia.org/wiki/Derangement
%
%   See also RANDPERM, RAND, NCHOOSEK, PERMS
%            COMBN, SHAKE, RANDSWAP (Matlab File Exchange)
% version 4.0 (may 2019)
% tested in Matlab 2018a
% (c) Jos van der Geest
% Matlab File Exchange Author ID: 10584
% email: samelinoa@gmail.com
% File history:
% 1.0 (nov 2005) - created
% 2.0 (mar 2009) - made suitable for the File Exchange
% 2.1 (mar 2009) - a few textual changes
% 2.2 (dec 2010) - randpermfull did not return all possible permutations,
%  only the cyclic (2) permutations (credits to Jan Simon & Derek O'Connor).
%  A rejection algorithm generally has a good performance, although it is
%  less elegant than the (slow) algorithm proposed by Martinez (SIAM 2008).
% 3.0 (feb 2016) - revived, implemented a for-loop to check for a proper
%  derangement, rather then using ANY.
% 4.0 (may 2019) - fixed serious error, where last element could be the
% same as its index (Thanks to Vincent L for pointing this out!)
% % a minimal error check for common mistakes
if nargin ~= 1 || numel(N) ~= 1 || ~isnumeric(N) || fix(N) ~= N || N < 2
    error('Argument should be a positive integer larger than 1') ;
end
if N == 2
    % the trivial case of two elements
    R = [2 1] ;
else
    % a fast rejection scheme
    hasViolation = true ; 
    while hasViolation
        R = randperm(N) ;      % create a possible derangement
        hasViolation = false ; % and assume that it is :-)
        for k = 1:N % check for derangement by looping over the elements
            if R(k) == k % violation found
                hasViolation = true ; % try again
                break ; % no need to check further
            end
        end
    end
end
% Old cyclic algorithm, not producing all possible derangements
% K = [1:N] + ceil((N - R) .* rand(1,N))  ;
% % K(j) is a value on the interval [j+1, N] (j=1:N-1), K(end) is N+1
% for j=1:N-1, "R(j) <-> R(K(j))"
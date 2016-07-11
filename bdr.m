% function [stats] = bdr(X,Y,A,B)
% Bidimensional Regression: A Matlab Implementation by Thomas J. Pingel
%
% This matlab function is a coding of Friedman and Kohler's 
% (2003) Euclidean bidimensional regression solution.
%
% References: 
% Friedman, A., & Kohler, B. (2003). Bidimensional Regression: Assessing
%   the Configural Similarity and Accuracy of Cognitive Maps and Other Two-Dimensional 
%   Data Sets. Psychological Methods, 8(4), 468-491.
%   Available on the web: 
%   http://www.psych.ualberta.ca/~alinda/PDFs/Friedman%20Kohler%20%5B03-Psych%20Methods%5D.pdf
% Nakaya, T. (1997). Statistical inferences in bidimensional regression
%   models. Geographical Analysis, 29, 169-186.
% Tobler, W. (1994). Bidimensional Regression. Geographical Analysis, 26,
%   187-212.
%
% This code copyrighted March 19, 2008.
%
% Free for use, but please cite if used for data prepared for publication.
%
% Thomas. J Pingel
% Department of Geography
% 1832 Ellison Hall 
% University of California Santa Barbara
% Santa Barbara, CA 93106-4060
% pingel@geog.ucsb.edu
%
% The function can be called with zero, one, or four arguments.
% Zero - Uses sample a dataset provided in article - e.g., [stats] = bdm();
% One  - User data in matrix form - e.g., [stats] = bdm([X Y A B]);
% Four - User data in (n-by-1) vector form - e.g., [stats] = bdm(X,Y,A,B);
%
% Any rows with any NaN will be eliminated.  If the missing values must be estimated
% beforehand, this must be done prior to passing them to the bdr function.
%
% The output is a structure containing the parameters of the regression and distortion.  These are:
%
% n -           Number of points used in the analysis 
% beta1 -       Regression parameter
% beta2 -       Regression parameter
% alpha1 -      Translation dispacement for X
% alpha2 -      Translation displacement for Y
% scale (phi) - A value greater than one means dependent is larger (expansion), 
%               a value smaller than one indicates dependent coordinates are 
%               smaller (contraction).
% theta -       A measure of angular displacement.  Negative indicates a
%               clockwise rotation is necessary.
% aPrime -      Predicted X coordinates
% bPrime -      Predicted Y coordinates
% rsquare -     Overall amount of variance explained by regression
% D -           Square root of the unexplained variance of between
%               dependent and predicted values (AB and ABprime)
% Dmax -        Maximum value of D.
% DI -          Distortion Index.  A measure of unexplained variance.
% F -           F-ratio, from Nakaya, 1997, Equation 50
% P -           p value for F, also from Nakaya.
% 
% Of these, alphas, phi, and theta are likely of most use to the researcher
% as specific descriptions of parts of the distortion.
% Rsquare and DI are likely the most useful as overall or general measures of
% distortion.
%
% Updates:
% March 19 - Added a step to remove any rows containing NaNs

function [stats] = bdr(varargin)
if length(varargin)==0 | length(varargin)==1 | length(varargin)==4
    if length(varargin)==0
        disp('Bidimensional Regression Function');
        disp('Running Case 1 test sample published in Friedman and Kohler, 2003');
        matrix = [0 0 12 16; 10 0 19 19; 8 5 21 18; 6 5 18 14]; % Last row to test elimination of NaNs
        disp('     X     Y     A     B');
        disp(matrix);
        X = matrix(:,1);
        Y = matrix(:,2);
        A = matrix(:,3);
        B = matrix(:,4);
    elseif length(varargin)==1
        % Assume matrix input
        matrix = varargin{1,1};
        X = matrix(:,1);
        Y = matrix(:,2);
        A = matrix(:,3);
        B = matrix(:,4);
    elseif length(varargin)==4
        X = varargin{1,1};
        Y = varargin{1,2};
        A = varargin{1,3};
        B = varargin{1,4};
        if (size(X,2) > size(X,1)) % Attempts a quick fix if vectors are in (1xN) instead of (Nx1)
            X=X';
            Y=Y';
            A=A';
            B=B';
        end
    end
    
    % Eliminate rows with NaNs
    badRows = find(isnan(X) | isnan(Y) | isnan(A) | isnan(B));
    X(badRows)=[];
    Y(badRows)=[];
    A(badRows)=[];
    B(badRows)=[];

    % These equations are given in Table 2, page 475 of the article.
    beta1 = (sum((X-mean(X)).*(A-mean(A))) + sum((Y-mean(Y)).*(B-mean(B))))/(ssq(X)+ssq(Y));
    beta2 = (sum((X-mean(X)).*(B-mean(B))) - sum((Y-mean(Y)).*(A-mean(A))))/(ssq(X)+ssq(Y));
    scale = (beta1^2 + beta2^2)^.5;
    theta = 180*atan2(beta2,beta1)/pi;
    alpha1 = mean(A) - beta1*mean(X) + beta2*mean(Y);
    alpha2 = mean(B) - beta2*mean(X) - beta1*mean(Y);
    aPrime = alpha1 + beta1*X - beta2*Y;
    bPrime = alpha2 + beta2*X + beta1*Y;
    rsquare = 1 - sum((A-aPrime).^2 + (B-bPrime).^2)/sum(ssq(A)+ssq(B));
    D = sqrt(sum((A-aPrime).^2 + (B-bPrime).^2));    % Equivalent to D_AB in the article
    Dmax = sqrt(ssq(A) + ssq(B));                % Equivalent to Dmax_AB
    % Dmax = sqrt(ssq(X) + ssq(Y));                % Note that these are square roots of equivalent numbers in table 3 of paper
    DI = sqrt(1-rsquare);                        % This is how DI is computed in sample and equation 12, page 479
    F = (2*size(X,1)-4) * rsquare / (2*(1-rsquare));
    P = 1-(fcdf(F,2,2*(size(X,1))-4));

    % Assemble results into a structure.  Sure, I could have done it
    % earlier, without using other variables.  This way I can easily
    % control the final order they appear in.
    stats.n = length(A); % Number of rows in analysis
    stats.beta1 = beta1;
    stats.beta2 = beta2;
    stats.alpha1 = alpha1;
    stats.alpha2 = alpha2;
    stats.scale = scale;
    stats.theta = theta;
    stats.aPrime = aPrime;
    stats.bPrime = bPrime;
    stats.rsquare = rsquare;
    stats.D = D;
    stats.Dmax = Dmax;
    stats.DI = DI;
    stats.F = F;
    stats.P = P;
    
    
else % Give instructions if the program wasn't called with the right numbers of parameters.
    disp('Insufficient number of parameters.');
    disp('Calling function with no parameters runs test case.');
    disp('Calling function with one parameter assumes matrix entry - [X,Y,A,B]');
    disp('Calling function with four parameters assumes seperate entry - X,Y,A,B');
end            



function [output] = ssq(input)
   output=sum((input-mean(input)).^2);
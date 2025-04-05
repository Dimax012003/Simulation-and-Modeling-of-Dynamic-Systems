function [theta,y_hat]=least_squares(X,y)
    theta=inv(X'*X)*X'*y;
    y_hat=X*theta;
end
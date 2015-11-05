function [cost]= GenCost(x,A,B,C)

cost=A.*x^2+B.*x+C;
end
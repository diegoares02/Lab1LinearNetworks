function [epsilon,epsilon_r]=get_epsilon(R,N,Pw,Fw,RL)
if length(R)==N
    X=(1/sqrt((10^(RL/10))-1))*abs(sum(Pw)/sum(Fw));
    epsilon=sqrt(X^2+1);
    epsilon_r=epsilon/sqrt(epsilon^2-1);
else
    epsilon_r=1;
    epsilon=(1/(sqrt(10^(RL/10)-1)))*(sum(Pw)/sum(Fw));
end
fprintf('Epsilon: %f Epsilon R: %f\n',epsilon,epsilon_r);
end
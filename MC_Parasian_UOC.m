function [ Price ] = MC_Parasian_UOC( S, K, t, T, n, sigma, r, N, B, tau)
%Computation of the MonteCarlo price of a Up and Out Parasian Call option
%
%Input:
%S Underlying stock price
%K Strike price
%t Actual time
%T Maturity time 
%n number of equally spaced time intervals
%sigma Standard deviation
%r Risk free rate
%N number of simulations
%B Barrier
%
%Output:
%Price

S(1)=S; %Starting point of GBM
dt=(T-t)/n; %Time interval

Sample_Payoff=zeros(1,N); %Preallocation for efficiency
Payoff=zeros(1,N); %Preallocation for efficiency

counting_time_over_barrier=zeros(1,N); %Counting variable for the time over the barrier
root_dt=sqrt(dt);


%Sample paths generator via Euler SDE approximation

for j=1:N
    for i=2:n+1
        S(i)=S(i-1)+r*S(i-1)*dt+sigma*S(i-1)*root_dt*randn;  
        if S(i)>=B
            counting_time_over_barrier(j)=counting_time_over_barrier(j)+dt;
        end
    end
    Sample_Payoff(j)=S(n+1);
end

%Up and Out condition

for j=1:N 
    if counting_time_over_barrier(j)>=tau
        Payoff(j)=0;
    else
        Payoff(j)=max(Sample_Payoff(j)-K,0);
    end
end

Price=exp(-r*(T-t))*mean(Payoff);

end
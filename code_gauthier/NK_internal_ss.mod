close all;

@#define ramsey = 1

var pi, c, h, r, a, y;
varexo u;
parameters beta, rho, epsilon, kappa_P, chi, gamma, H;

beta	= 0.995;
gamma	= 3;
kappa_P	= 50;
epsilon	= 8;
chi		= 1;
rho		= 0.95;
H		= 1;
chi 	= ((epsilon-1))/epsilon/(H*H^(gamma));


model;
	@#if ramsey == 0
		r-steady_state(r)=1.5*(pi-steady_state(pi));
	@#endif

	1/c = beta*r/(c(+1)*pi(+1));
	kappa_P*pi*(pi-1) = (1-epsilon) + epsilon*c*chi*h^(gamma)/exp(a) +  beta*c/c(+1)*y(+1)/y * kappa_P*( pi(+1) - 1 )*pi(+1);
	y = c+(kappa_P/2)*(pi-1)^2*y;
	y = exp(a)*h;
	a = rho*a(-1)+u;
end;

@#if ramsey == 1
	initval;
		r = 1/beta+0.01;
	end;
@#endif




steady_state_model;
	% write the steady state with r as given
	@#if ramsey == 0
		r = 1/beta+0.01;
	@#endif
	pi = beta*r;
	h = ((kappa_P*pi*(pi-1)*(1-beta)-(1-epsilon))/((1-(kappa_P/2)*(pi-1)^2)*epsilon*chi))^(1/(1+gamma));
	y = h;
	c = (1-(kappa_P/2)*(pi-1)^2)*y;
	a = 0;
end;

shocks;
 var u;stderr 0.01;
end;

@#if ramsey == 1
	planner_objective(log(c)-chi*((h^(1+gamma))/(1+gamma)));
	ramsey_model(planner_discount = 0.95,instruments=(r));
@#endif

%options_.ramsey_policy=1;
steady;
resid;
stoch_simul(order=1) y pi r;

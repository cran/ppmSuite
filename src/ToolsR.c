#include <R.h>
#include <Rmath.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////
//Random Numbers
///////////////////////////////////////////////////////////////////////////////////////////////////////

void F77_SUB(rndstart)(void) { GetRNGstate(); }

void F77_SUB(rndend)(void) { PutRNGstate(); }

double F77_SUB(betar)(double *shape1, double *shape2)
{
	return rbeta(*shape1, *shape2);
}

double F77_SUB(chisqr)(double *df)
{
	return rchisq(*df);
}

double F77_SUB(gammar)(double *alpha, double *beta)
{
	return rgamma(*alpha, 1.0)/(*beta);
}

double F77_SUB(normr)(double *mean, double *sd)
{
	return rnorm(*mean, *sd);
}

double F77_SUB(unifr)(double *a, double *b)
{
	return runif(*a, *b);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Probability Density Functions
///////////////////////////////////////////////////////////////////////////////////////////////////////

double F77_SUB(betad)(double *x, double *shape1, double *shape2, int *give_log)
{
	return dbeta(*x, *shape1, *shape2, *give_log);
}

double F77_SUB(binomd)(double *x, double *size, double *prob, int *give_log)
{
	return dbinom(*x, *size, *prob, *give_log);
}

double F77_SUB(gammad)(double *x, double *shape, double *scale, int *give_log)
{
	return dgamma(*x, *shape, *scale, *give_log);
}

double F77_SUB(lnormd)(double *x, double *meanlog, double *sdlog, int *give_log)
{
	return dlnorm(*x, *meanlog, *sdlog, *give_log);
}

double F77_SUB(normd)(double *x, double *mean, double *sd, int *give_log)
{
	return dnorm(*x, *mean, *sd, *give_log);
}

double F77_SUB(poisd)(double *x, double *lambda, int *give_log)
{
	return dpois(*x, *lambda, *give_log);
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
//Cumulative Distribution Functions
///////////////////////////////////////////////////////////////////////////////////////////////////////

double F77_SUB(betap)(double *x, double *shape1, double *shape2, int *lower_tail, int *log_p)
{
	return pbeta(*x, *shape1, *shape2, *lower_tail, *log_p);
}

double F77_SUB(lnormp)(double *x, double *meanlog, double *sdlog, int *lower_tail, int *log_p)
{
        return plnorm(*x, *meanlog, *sdlog, *lower_tail, *log_p);
}

double F77_SUB(normp)(double *x, double *mean, double *sd, int *lower_tail, int *log_p)
{
	return pnorm(*x, *mean, *sd, *lower_tail, *log_p);
}



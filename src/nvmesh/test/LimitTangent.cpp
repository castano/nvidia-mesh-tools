
#include <math.h>
#include <nvmath/nvmath.h>

#include <stdio.h>
/*
inline static float alphaWeightU(int i, int n)
{
	const float cp = cosf(PI/n);
	const float denom = (n * sqrtf(4 + cp*cp));
	return (1.0f / n + cp / denom) * cosf((2.0f * PI * i) / n);  
}

inline static float betaWeightU(int i, int n)
{
	const float cp = cosf(PI/n);
	const float denom = (n * sqrtf(4 + cp*cp));
	return (1.0f / denom) * cosf((2.0f * PI * i + PI) / n);  
}
*/
inline static float alphaWeightU(int i, int n)
{
	const float cp = cosf(PI/n);
	const float denom = (n * sqrtf(4 + cp*cp));
	return (1.0f / n + cp / denom) * cosf((2.0f * PI * i) / n);  
}

inline static float betaWeightU(int i, int n)
{
	const float cp = cosf(PI/n);
	const float denom = (n * sqrtf(4 + cp*cp));
	return (1.0f / denom) * cosf((2.0f * PI * i + PI) / n);  
}

inline static float alphaWeightV(int i, int n)
{
	const float cp = cosf(PI/n);
	const float denom = (n * sqrtf(4 + cp*cp));
	return (1.0f / n + cp / denom) * sinf((2.0f * PI * i) / n);  
}

inline static float betaWeightV(int i, int n)
{
	const float cp = cosf(PI/n);
	const float denom = (n * sqrtf(4 + cp*cp));
	return (1.0f / denom) * sinf((2.0f * PI * i + PI) / n);  
}



int main(void)
{
	float total = 0;
	const int n = 3;
	for(int i = 0; i < n; i++)
	{
		float alpha = alphaWeightU(i, n);
		float beta = betaWeightU(i, n);
		total += alpha + beta;
		printf("alpha %i = %f\n", i, alpha);
		printf("  beta %i = %f\n", i, beta);
	}
	printf("\ntotal = %f\n\n", total);

	total = 0;
	for(int i = 0; i < n; i++)
	{
		float alpha = alphaWeightV(i, n);
		float beta = betaWeightV(i, n);
		total += alpha + beta;
		printf("alpha %i = %f\n", i, alpha);
		printf("  beta %i = %f\n", i, beta);
	}
	printf("\ntotal = %f\n\n", total);

	return 0;
}


#include "device_lib.cuh"


__device__  void CreateQ_By_float(struct_Q_T *out, float x, float y, float z, float w)
{
	out->x = x;
	out->y = y;
	out->z = z;
	out->w = w;
}

__device__  float  Get_QNorm(struct_Q_T *Q)
{
	return sqrtf(Q->x*Q->x + Q->y*Q->y + Q->z*Q->z + Q->w*Q->w);
}

__device__ void Get_QPow(struct_Q_T *Q, float pow)
{
	float A = Get_QNorm(Q);
	float theta = 0.0f;
	float B = 0.0f;
	float R = 0.0f;
	if (pow > 0.0f && A>0.000001f)
	{
		float coef = 1.0f;
		if (A<1.0f)
		{
			//printf("%f *******\n", A);
			coef = 1 / A;
			Q->x *= coef;
			Q->y *= coef;
			Q->z *= coef;
			Q->w *= coef;

		}
		A = Get_QNorm(Q);
		//printf("%f +++++++++\n", A);
		theta = acosf(Q->w / A)*pow;
		B = sqrt(A*A - Q->w*Q->w);
		R = exp2f(logf(A / coef)* pow);
		Q->x = R*sinf(theta)*(Q->x / B);
		Q->y = R*sinf(theta)*(Q->y / B);
		Q->z = R*sinf(theta)*(Q->z / B);
		Q->w = R*cosf(theta);

	}
	else
	{
		//printf("%f --------\n", A);
		Q->w = 0.0f;
		Q->x = 0.0f;
		Q->y = 0.0f;
		Q->z = 0.0f;

	}
}
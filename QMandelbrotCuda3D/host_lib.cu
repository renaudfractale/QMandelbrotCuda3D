#include "host_lib.cuh"


__host__  void CreateQ_By_float_H(struct_Q_T *out, float x, float y, float z, float w)
{
	out->x = x;
	out->y = y;
	out->z = z;
	out->w = w;
}

__host__  float  Get_QNorm_H(struct_Q_T *Q)
{
	return sqrtf(Q->x*Q->x + Q->y*Q->y + Q->z*Q->z + Q->w*Q->w);
}

__host__ void Get_QPow_H(struct_Q_T *Q, float pow)
{
	float A = Get_QNorm_H(Q);
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
		A = Get_QNorm_H(Q);
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

__host__ int  GetQIter_H(struct_P_Simulation_T *P_Simulation_DEVICE, int  *x_filter, int  *y_filter, int *z_filter, int *w_filter)
{
	//int Tempindex = 0;
	struct_Q_T Q_Current;
	float w, x, y, z;
	int iter_computed;
	//X
	x = ((float)*x_filter)*P_Simulation_DEVICE->X.step + P_Simulation_DEVICE->X.start;
	//Y
	y = ((float)*y_filter)*P_Simulation_DEVICE->Y.step + P_Simulation_DEVICE->Y.start;
	//Z
	z = ((float)*z_filter)*P_Simulation_DEVICE->Z.step + P_Simulation_DEVICE->Z.start;
	//W
	w = ((float)*w_filter)*P_Simulation_DEVICE->W.step + P_Simulation_DEVICE->W.start;



	CreateQ_By_float_H(&Q_Current, x, y, z, w);

	for (iter_computed = 0; iter_computed <= P_Simulation_DEVICE->Iter.max; iter_computed++)
	{
		Get_QPow_H(&Q_Current, P_Simulation_DEVICE->power);
		Q_Current.x += x;
		Q_Current.y += y;
		Q_Current.z += z;
		Q_Current.w += w;

		if (Get_QNorm_H(&Q_Current) > P_Simulation_DEVICE->rMax)
		{
			if (iter_computed > 0)
				iter_computed--;
			return iter_computed;
		}
	}
	if (iter_computed > 0)
		iter_computed--;
	return iter_computed;
}
__host__ bool  FilterQ_H(int *Filter, int *Nx, int *Ny, int *Nz, int *Nw, int iter,	struct_P_Simulation_T *P_Simulation)
{
	if (*Filter == 0)
		return true;
	int iter_computed = 0;

	int pasx = 1;
	if (P_Simulation->X.NbStep == 1)
		pasx = 0;

	int pasy = 1;
	if (P_Simulation->Y.NbStep == 1)
		pasy = 0;

	int pasz = 1;
	if (P_Simulation->Z.NbStep == 1)
		pasz = 0;

	int pasw = 1;
	if (P_Simulation->W.NbStep == 1)
		pasw = 0;


	for (int x_filter = *Nx - pasx; x_filter <= *Nx + pasx; x_filter++)
	{
		for (int y_filter = *Ny - pasy; y_filter <= *Ny + pasy; y_filter++)
		{
			for (int z_filter = *Nz - pasz; z_filter <= *Nz + pasz; z_filter++)
			{
				for (int w_filter = *Nw - pasw; w_filter <= *Nw + pasw; w_filter++)
				{
					iter_computed = GetQIter_H(P_Simulation, &x_filter, &y_filter, &z_filter, &w_filter);
					if (*Filter == 1)
					{
						if (iter_computed != iter)
							return true;
					}
					else //filter==2
					{
						if (iter_computed == 0)
							return true;
					}
				}
			}
		}
	}
	return false;
}

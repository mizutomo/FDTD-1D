#ifndef __MUR_1ST_H__
#define __MUR_1ST_H__

void ex_boundary(double** ex, double** hist_ex, double* dy, 
								 int nx, int ny, double dt);
void ey_boundary(double** ey, double** hist_ey, double* dx,
								 int nx, int ny, double dt);

#endif //__MUR_1ST_H__

// Mur1次吸収境界
#include <stdio.h>
#include <stdlib.h>
#include "common.h"

void hist_update_x(double** ex, double** hist_ex, int nx, int ny);
void hist_update_y(double** ey, double** hist_ey, int nx, int ny);

void ex_boundary(double** ex, double** hist_ex, double* dy, 
								 int nx, int ny, double dt)
{
	int i;
	double c = LIGHT * dt;

	for (i = 0; i < nx; i++) {
		ex[i][0] = hist_ex[i][1] 
			+ (c - dy[0])    / (c + dy[0])    * (ex[i][1]    - hist_ex[i][0]);
		ex[i][ny] = hist_ex[i][2]
			+ (c - dy[ny-1]) / (c + dy[ny-1]) * (ex[i][ny-1] - hist_ex[i][3]);
	}
}

void ey_boundary(double** ey, double** hist_ey, double* dx,
								 int nx, int ny, double dt)
{
	int j;
	double c = LIGHT * dt;

	for (j = 0; j < ny; j++) {
		ey[0][j]  = hist_ey[j][1]
			+ (c - dx[0]) / (c + dx[0]) * (ey[1][j] - hist_ey[j][0]);
		ey[nx][j] = hist_ey[j][2]
			+ (c - dx[nx-1]) / (c + dx[nx-1]) * (ey[nx-1][j] - hist_ey[j][3]);
	}
}

void hist_update(double** ex, double** ey,
								 double** hist_ex, double** hist_ey,
								 int nx, int ny)
{
	hist_update_x(ex, hist_ex, nx, ny);
	hist_update_y(ey, hist_ey, nx, ny);
}

void hist_update_x(double** ex, double** hist_ex, int nx, int ny)
{
	int i;

	for (i = 0; i < nx; i++) {
		hist_ex[i][0] = ex[i][0];
		hist_ex[i][1] = ex[i][1];
		hist_ex[i][2] = ex[i][ny-1];
		hist_ex[i][3] = ex[i][ny];
	}
}

void hist_update_y(double** ey, double** hist_ey, int nx, int ny)
{
	int j;

	for (j = 0; j < ny; j++) {
		hist_ey[j][0] = ey[0][j];
		hist_ey[j][1] = ey[1][j];
		hist_ey[j][2] = ey[nx-1][j];
		hist_ey[j][3] = ey[nx][j];
	}
}

// FDTDメインプログラム
// 10.01.06 by T.Mizukusa
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "normal.h"
//#include "ade.h"
//#include "mur_1st.h"
#include "berengerPML.h"
//#include "fdtd.h"

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define _(X) (X) && (X)

// アルゴリズムを現すEnumerate Type
enum eALG {
	NORMAL,   // Normal-FDTD
	ADE,      // ADE(平均化する)-FDTD
	ADE2      // ADE(平均化しない)-FDTD
};

// グローバル変数
// 材料定数の宣言
mat_t material[2] = {
	{1.0, 1.0, 0.0},      // 真空
	{1.0, 1.0, 100e+6},   // PEC
};

// X, Y座標
typedef struct _point_t {
	double x;
	double y;
} point_t;

// ベクトルデータ型
typedef struct _vec_t {
	double xvec;
	double yvec;
} vec_t;

// 長方形データ型
typedef struct _rectangle_t {
	point_t points[4];
} rectangle_t;

int total_nx, total_ny;

double *dx, *dy;

conform_t** mat;

double **ex, **ey;
double **hz;
double *hzx, *hzy;

double **cex, **cexly;
double **cey, **ceylx;
double **chzlx, **chzly;

double dt;

int bck;
int bck_calib;

// プロトタイプ宣言
double** alloc_2d_array_dbl(int nx, int ny);
void initialize_2d_array_dbl(double** ary, int nx, int ny);
conform_t** alloc_2d_array_conform_t(int nx, int ny);
void initialize_2d_array_conform_t(conform_t** ary, int nx, int ny);
double calc_stimulus(int step, double dt);
void print_point_value_header(FILE* fp);
void print_point_value(FILE* fp, double** data, double time, double stimulus);
void print_line_value(FILE* fp, double** data, int step, int mode);
void print_map_value(FILE* fp, double** data, int step);
void set_cylinder_area_ratio(int i, int j, double x, double y, double dx, double dy, double cx, double cy, double r);
void set_cylinder_lx_ratio(int i, int j, double x, double y, double dx, double dy, double cx, double cy, double r);
void set_cylinder_ly_ratio(int i, int j, double x, double y, double dx, double dy, double cx, double cy, double r);
float get_cylinder_xsection_area(double x, double y, double dx, double dy, double cx, double cy, double r);
void set_rectangle_area_ratio(int i, int j, double x, double y, double dx, double dy, rectangle_t rectangle);
float get_rectangle_xsection_area(double x, double y, double dx, double dy, rectangle_t rectangle);
void set_rectangle_lx_ratio(int i, int j, double x, double y, double dx, double dy, rectangle_t rectangle);
void set_rectangle_ly_ratio(int i, int j, double x, double y, double dx, double dy, rectangle_t rectangle);
point_t get_cross_point_rectangle(double x0, double y0, double x1, double y1, rectangle_t rectangle, int dir);
double get_value_from_linear_equation(double x, point_t p0, point_t p1);
float find_max_edge_ratio(int i, int j);
float find_permittivity_ratio(int i, int j);
void rotate_point(float radian, point_t* point);
int judge_inside_rectangle(double x, double y, rectangle_t rectangle);
int get_cross_point_with_two_line(point_t grid0, point_t grid1, point_t seg0, point_t seg1, 
																	int dir, point_t* xpoint);



// 二乗関数
double square(double x)
{
	return x * x;
}

// 最小のセルを見つける関数
double find_min_cell(double* ary, int leng)
{
	int i;
	double min;

	min = 1e+30;
	for (i = 0; i < leng; i++) {
		if (ary[i] <= min) {
			min = ary[i];
		}
	}

	return min;
}

// 配列領域の確保&初期化
void alloc_memory_for_analysis_area()
{
	ex = alloc_2d_array_dbl(    total_nx, total_ny+1);
	initialize_2d_array_dbl(ex, total_nx, total_ny+1);

	ey = alloc_2d_array_dbl(    total_nx+1, total_ny);
	initialize_2d_array_dbl(ey, total_nx+1, total_ny);

	hz = alloc_2d_array_dbl(    total_nx, total_ny);
	initialize_2d_array_dbl(hz, total_nx, total_ny);

	cex = alloc_2d_array_dbl(    total_nx, total_ny+1);
	initialize_2d_array_dbl(cex, total_nx, total_ny+1);

	cexly = alloc_2d_array_dbl(    total_nx, total_ny+1);
	initialize_2d_array_dbl(cexly, total_nx, total_ny+1);

	cey = alloc_2d_array_dbl(    total_nx+1, total_ny);
	initialize_2d_array_dbl(cey, total_nx+1, total_ny);

	ceylx = alloc_2d_array_dbl(    total_nx+1, total_ny);
	initialize_2d_array_dbl(ceylx, total_nx+1, total_ny);

	chzlx = alloc_2d_array_dbl(    total_nx, total_ny);
	initialize_2d_array_dbl(chzlx, total_nx, total_ny);

	chzly = alloc_2d_array_dbl(    total_nx, total_ny);
	initialize_2d_array_dbl(chzly, total_nx, total_ny);
}

// X,Y平面に対して、グリッド間隔を設定(X軸)
void set_delta_x()
{
	int i;
	for (i = 0; i < total_nx; i++) dx[i] = DX;
}

// X,Y平面に対して、グリッド間隔を設定(Y軸)
void set_delta_y()
{
	int i;
	for (i = 0; i < total_ny; i++) dy[i] = DY;
}

// CFL値の算出
void set_cfl_constant(float cfl, double* min_dx, double* min_dy)
{
	(*min_dx) = find_min_cell(dx, total_nx);
	(*min_dy) = find_min_cell(dy, total_ny);

	dt = cfl * 1.0 / 
		(LIGHT * sqrt(square(1.0/(*min_dx)) + 
									square(1.0/(*min_dy))));
}

// BCKセルの面積比をタイムステップに応じて修正
void calibrate_bck_cell_area(float cfl_bck)
{
	int i, j;
	float amin_ratio;

	for (i = 0; i < NX; i++) {
		for (j = 0; j < NY; j++) {
			if (mat[i][j].conform_flag == 1) {
				// (1) Initialize
				amin_ratio = square(cfl_bck / 1.0);

				// (2) check the edge ratios
				amin_ratio = amin_ratio * find_max_edge_ratio(i, j);

				// (3) check permittivity
				amin_ratio = amin_ratio * find_permittivity_ratio(i, j);

				// (4) check cellsize
				// 不明確なので未実装

				// (5) Aratioの修正
				if (mat[i][j].a_ratio <= amin_ratio) {
					printf("Modify Aratio From %f to %f at [%d, %d]\n", mat[i][j].a_ratio, amin_ratio, i, j);
					mat[i][j].a_ratio = amin_ratio;
				}
			}
		}
	}
}

// 各セルの最大辺長比を見つける
float find_max_edge_ratio(int i, int j) 
{
	return MAX(mat[i][j].lx_ratio, mat[i][j].ly_ratio);
}

// 各セルのセカンドマテリアルの比誘電率を見つける
float find_permittivity_ratio(int i, int j)
{
	return material[mat[i][j].second_material].eps;
}

// XY2重配列のメモリ確保
// ヒープに領域を作りたいだけなので、freeはしてない。
double** alloc_2d_array_dbl(int nx, int ny)
{
	double** ary;
	int i;

	ary = (double**)malloc(sizeof(double*) * nx);
	for(i = 0; i < nx; i++) {
		ary[i] = (double*)malloc(sizeof(double) * ny);
	}

	return ary;
}

// XY2重配列のメモリ確保(CONFORM_T)
// ヒープに領域を作りたいだけなので、freeはしてない。
conform_t** alloc_2d_array_conform_t(int nx, int ny)
{
	conform_t** ary;
	int i;

	ary = (conform_t**)malloc(sizeof(conform_t*) * nx);
	for(i = 0; i < nx; i++) {
		ary[i] = (conform_t*)malloc(sizeof(conform_t) * ny);
	}

	return ary;
}

// 解析領域の斜め線フィル
void set_rectangle_material()
{
	int i, j;
	double radian;
	rectangle_t rectangle;
	double x, y;

	// 回転前の長方形の座標
	rectangle.points[0].x = 0.01; rectangle.points[0].y = 0.055;
	rectangle.points[1].x = 0.01; rectangle.points[1].y = 0.045;
	rectangle.points[2].x = 0.09; rectangle.points[2].y = 0.045;
	rectangle.points[3].x = 0.09; rectangle.points[3].y = 0.055;

	// 長方形を45度(PI/4)回転
	radian = 2*PI/360 * 80;
	for (i = 0; i < 4; i++) {
		rotate_point(radian, &(rectangle.points[i]));
		printf("%f, %f, 0.0\n", rectangle.points[i].x, rectangle.points[i].y);
	}

	//start_x = start_y = 0.0;
	//for (i = 0; i < L; i++) start_x += dx[i];
	//for (j = 0; j < L; j++) start_y += dy[j];

	x = 0.0;
	for (i = 0; i < NX; i++) {
		y = 0.0;
		for (j = 0; j < NY; j++) {
			set_rectangle_area_ratio(i, j, x, y, dx[i+L], dy[j+L], rectangle);
			set_rectangle_lx_ratio(i, j, x, y, dx[i+L], dy[j+L], rectangle);
			set_rectangle_ly_ratio(i, j, x, y, dx[i+L], dy[j+L], rectangle);
			y += dy[j+L];
		}
		x += dx[i+L];
	}
}

// 微小領域と長方形の重なり具合を算出
void set_rectangle_area_ratio(int i, int j, double x, double y, double dx, double dy, rectangle_t rectangle)
{
	int f0, fx, fy, fxy; 
	double ratio;

	f0  = judge_inside_rectangle(x,    y,    rectangle);
	fx  = judge_inside_rectangle(x+dx, y,    rectangle);
	fy  = judge_inside_rectangle(x,    y+dy, rectangle);
	fxy = judge_inside_rectangle(x+dx, y+dy, rectangle);

	if (f0 == 1 && fx == 1 && fy == 1 && fxy == 1) {
		// グリッド全ての点が、長方形内に存在する。
		mat[i][j].master_material = 1;
		mat[i][j].a_ratio = 0.0;
		mat[i][j].conform_flag = 0;
	} else if (f0 == 0 && fx == 0 && fy == 0 && fxy == 0) {
		// グリッド全ての点が、長方形外に存在する。
		mat[i][j].master_material = 0;
		mat[i][j].a_ratio = 1.0;
		mat[i][j].conform_flag = 0;
	} else {
		ratio = get_rectangle_xsection_area(x, y, dx, dy, rectangle);
		mat[i][j].a_ratio = 1.0 - ratio;
		if (bck == 0) {    // BCKが選択されなかった場合、階段近似する。
			if (ratio >= 0.5) {
				mat[i][j].master_material = 1;   // 1グリッド内で50%以上導体が存在する時、マスターの材質をPECに設定
			} else {
				mat[i][j].master_material = 0;   // 1グリッド内で導体の割合が50%以下の場合、マスターの材質を真空に設定
			}
			mat[i][j].conform_flag = 0;
		} else {           // BCKセルを選択した場合
			mat[i][j].master_material = 1;
			mat[i][j].second_material = 0;
			mat[i][j].conform_flag = 1;
		}
	}
}

// X軸と長方形との交点を解析的に算出
void set_rectangle_lx_ratio(int i, int j, double x, double y, double dx, double dy, rectangle_t rectangle)
{
	int f1, f2;    // 長方形の内外フラグ
	point_t xpoint;

	f1 = judge_inside_rectangle(x,    y, rectangle);
	f2 = judge_inside_rectangle(x+dx, y, rectangle);

	if (f1 == 1 && f2 == 1) {
		// dxの線分が完全に長方形の中に存在する場合
		mat[i][j].lx_ratio = 0.0;
	} else if (f1 == 0 && f2 == 0) {
		// dxの線分が完全に長方形の外に存在する場合
		mat[i][j].lx_ratio = 1.0;
	} else if (f1 == 0 && f2 == 1) {
		xpoint = get_cross_point_rectangle(x, y, x+dx, y+dy, rectangle, 0);
		mat[i][j].lx_ratio = (xpoint.x - x) / dx;
	} else if (f1 == 1 && f2 == 0) {
		xpoint = get_cross_point_rectangle(x, y, x+dx, y+dy, rectangle, 0);
		mat[i][j].lx_ratio = (x + dx - xpoint.x) / dx;
	} else {
		fprintf(stderr, "[Error] Invalid cross seciton: X-Axis.\n");
		exit(1);
	}
}

// Y軸と長方形との交点を解析的に算出
void set_rectangle_ly_ratio(int i, int j, double x, double y, double dx, double dy, rectangle_t rectangle)
{
	int f1, f2;    // 長方形の内外フラグ
	point_t xpoint;

	f1 = judge_inside_rectangle(x, y,    rectangle);
	f2 = judge_inside_rectangle(x, y+dy, rectangle);

	if (f1 == 1 && f2 == 1) {
		// dyの線分が完全に長方形の中に存在する場合
		mat[i][j].ly_ratio = 0.0;
	} else if (f1 == 0 && f2 == 0) {
		// dyの線分が完全に長方形の外に存在する場合
		mat[i][j].ly_ratio = 1.0;
	} else if (f1 == 0 && f2 == 1) {
		xpoint = get_cross_point_rectangle(x, y, x+dx, y+dy, rectangle, 1);
		mat[i][j].ly_ratio = (xpoint.y - y) / dx;
	} else if (f1 == 1 && f2 == 0) {
		xpoint = get_cross_point_rectangle(x, y, x+dx, y+dy, rectangle, 1);
		mat[i][j].ly_ratio = (y + dy - xpoint.y) / dx;
	} else {
		fprintf(stderr, "[Error] Invalid cross seciton: Y-Axis.\n");
		exit(1);
	}
}

// 長方形の線分とグリッドの線分との交点を算出
point_t get_cross_point_rectangle(double x0, double y0, double x1, double y1, rectangle_t rectangle, int dir)
{
	int i, j;
	point_t xpoint, grid0, grid1;

	grid0.x = x0; grid0.y = y0;
	grid1.x = x1; grid1.y = y1;

	for (i = 0; i < 4; i++) {
		j = (i == 3)? 0 : i+1;
		if (get_cross_point_with_two_line(grid0, grid1, rectangle.points[i], rectangle.points[j], dir, &xpoint) 
				== 1) {
			return xpoint;
		}
	}

	return xpoint;
}

// 線分同士の交点を算出
int get_cross_point_with_two_line(point_t grid0, point_t grid1, point_t seg0, point_t seg1, 
																	int dir, point_t* xpoint)
{
	vec_t vec0, vec1;
	double lhs, rhs;
	double a;    // 1次関数の傾き(テンポラリ)
	point_t xp;  // 仮の交点
	double min, max;

	vec0.xvec = grid1.x - grid0.x; vec0.yvec = grid1.y - grid0.y;
	vec1.xvec = seg1.x  - seg0.x;  vec1.yvec = seg1.y  - seg0.y;

	// 平行性のチェック
	lhs = sqrt(square(vec0.xvec) + square(vec0.yvec)) * sqrt(square(vec1.xvec) + square(vec1.yvec));
	rhs = fabs(vec0.xvec * vec1.xvec + vec0.yvec * vec1.yvec);
	if (fabs(lhs - rhs) < 1e-12) {
		return 0;  // 二つのベクトルはほぼ平行
	}

	// 交点を算出
	a = (seg1.y - seg0.y) / (seg1.x - seg0.x);
	if (dir == 0) {
		// グリッド座標の最小値、最大値取得
		min = MIN(grid0.x, grid1.x);
		max = MAX(grid0.x, grid1.x);

		// X方向の交点を算出
		xp.y = grid0.y;
		xp.x = (xp.y - seg0.y) / a + seg0.x;

		// 交点がグリッド線分内にあるかチェック
		if (min < _(xp.x) < max) {
			xpoint->x = xp.x;
			xpoint->y = xp.y;
			return 1;
		} else {
			return 0;
		}
	} else {
		// グリッド座標の最小値、最大値取得
		min = MIN(grid0.y, grid1.y);
		max = MAX(grid0.y, grid1.y);

		// Y方向の交点を算出
		xp.x = grid0.x;
		xp.y = a * (xp.x - seg0.x) + seg0.y;

		// 交点がグリッド線分内にあるかチェック
		if (min < _(xp.y) < max) {
			xpoint->x = xp.x;
			xpoint->y = xp.y;
			return 1;
		} else {
			return 0;
		}

		return 1;
	}
}

// 微小領域と長方形との重なり具合をモンテカルロ法で算出
float get_rectangle_xsection_area(double x, double y, double dx, double dy, rectangle_t rectangle)
{
	int cycle, cycle_max, n;
	double mx, my;
	float area_ratio;

	n = 0;
	cycle_max = 1e+5;

	for (cycle = 0; cycle <= cycle_max; cycle++) {
		mx = x + (rand() / (RAND_MAX + 1.0)) * dx;
		my = y + (rand() / (RAND_MAX + 1.0)) * dy;

		if (judge_inside_rectangle(mx, my, rectangle) == 1) {
			n += 1;
		}
	}

	area_ratio = (float)n/(float)cycle_max;

	// 面積比が1 or 0を超えた場合に補正
	if (area_ratio > 1.0) area_ratio = 1.0;
	if (area_ratio < 0.0) area_ratio = 0.0;

	return area_ratio;
}

// ある点が長方形の中にあるか、外にあるかを判定
// 1: 中にある。0: 外にある。
int judge_inside_rectangle(double x, double y, rectangle_t rectangle)
{
	double y_prime[4];
	int i, j;

	for (i = 0; i < 4; i++) {
		j = (i == 3)? 0 : i + 1;
		y_prime[i] = get_value_from_linear_equation(x, rectangle.points[i], rectangle.points[j]);
	}

	if (y >= y_prime[0] && y >= y_prime[1] && y <= y_prime[2] && y <= y_prime[3]) {
		return 1;     // 点が長方形の中に入っている
	} else {
		return 0;
	}
}

// 1次関数からあるX座標のY座標値を取得
double get_value_from_linear_equation(double x, point_t p0, point_t p1)
{
	double x0, y0;
	double x1, y1;

	x0 = p0.x; y0 = p0.y;
	x1 = p1.x; y1 = p1.y;

	return (y1-y0) / (x1-x0) * (x-x0) + y0;
}

// 指定角度の回転
void rotate_point(float radian, point_t* point)
{
	point_t new_point;

	new_point.x = cos(radian) * (point->x - PX/2) - sin(radian) * (point->y - PY/2);
	new_point.y = sin(radian) * (point->x - PX/2) + cos(radian) * (point->y - PY/2);

	point->x = new_point.x + PX/2;
	point->y = new_point.y + PY/2;
}

// 解析領域の円柱フィル
void set_cylinder_material()
{
	int i, j;
	double x, y;
	double start_x, start_y;
	double cx, cy;
	double r = 10e-3;   // 10[mm]

	start_x = start_y = 0.0;
	for (i = 0; i < L; i++) start_x += dx[i];
	for (j = 0; j < L; j++) start_y += dy[j];

	cx = cy = 0.0;
	for (i = 0; i < (NX/2)+L; i++) cx += dx[i];
	for (j = 0; j < (NY/2)+L; j++) cy += dy[j];

	x = start_x;
	for (i = 0; i < NX; i++) {
		y = start_y;
		for (j = 0; j < NY; j++) {
			set_cylinder_area_ratio(i, j, x, y, dx[i+L], dy[j+L], cx, cy, r);
			set_cylinder_lx_ratio(i, j, x, y, dx[i+L], dy[j+L], cx, cy, r);
			set_cylinder_ly_ratio(i, j, x, y, dx[i+L], dy[j+L], cx, cy, r);
			y += dy[j+L];
		}
		x += dx[i+L];
	}
}

// 微小領域と円柱との重なり具合を算出
void set_cylinder_area_ratio(int i, int j, double x, double y, double dx, double dy, double cx, double cy, double r)
{
	double p0, px, py, pxy;
	double ratio;

	p0  = square(x    - cx) + square(y    - cy);
	px  = square(x+dx - cx) + square(y    - cy);
	py  = square(x    - cx) + square(y+dy - cy);
	pxy = square(x+dx - cx) + square(y+dy - cy);

	if (p0 <= square(r) && px <= square(r) && py <= square(r) && pxy <= square(r)) {
		// グリッドが全て導体に包まれている場合
		mat[i][j].master_material = 1;
		mat[i][j].a_ratio = 0.0;
		mat[i][j].conform_flag = 0;
	} else if (p0 > square(r) && px > square(r) && py > square(r) && pxy > square(r)) {
		// グリッドが全て誘電体の場合
		mat[i][j].master_material = 0;
		mat[i][j].a_ratio = 1.0;
		mat[i][j].conform_flag = 0;
	} else {
		ratio = get_cylinder_xsection_area(x, y, dx, dy, cx, cy, r);
		mat[i][j].a_ratio = 1.0 - ratio;
		if (bck == 0) {         // BCKが選択されなかった場合、階段近似する。
			if (ratio >= 0.5) {
				mat[i][j].master_material = 1;  // PECを設定
			} else {
				mat[i][j].master_material = 0;  // 真空を設定
			}
			mat[i][j].conform_flag = 0;
		} else {                // BCKセルを選択した場合
			mat[i][j].master_material = 1;    // 第1材料にPECを設定
			mat[i][j].second_material = 0;    // 第2材料に真空を設定
			mat[i][j].conform_flag = 1;
		}
	}
}

// 微小領域と円柱との重なり具合をモンテカルロ法で算出
float get_cylinder_xsection_area(double x, double y, double dx, double dy, double cx, double cy, double r)
{
	int cycle, cycle_max, n;
	double mx, my;
	float area_ratio;

	n = 0;
	//cycle_max = 1e+7;        //      乱数の発生回数
	cycle_max = 1e+5;        //      乱数の発生回数

	for(cycle = 0; cycle <= cycle_max; cycle++) {
		mx = x + (rand() /(RAND_MAX + 1.0)) * dx;
		my = y + (rand() /(RAND_MAX + 1.0)) * dy;

		if((square(mx - cx) + square(my - cy)) < square(r)) {
			n += 1;                 
		}
	}

	area_ratio = (float)n/(float)cycle_max;

	if (area_ratio > 1.0) area_ratio = 1.0;
	if (area_ratio < 0.0) area_ratio = 0.0;

	return area_ratio;
}

// X軸と円柱との交点を解析的に算出
void set_cylinder_lx_ratio(int i, int j, double x, double y, double dx, double dy, double cx, double cy, double r)
{
	double point1, point2;
	double xpoint;

	point1 = square(x - cx) + square(y - cy);
	point2 = square(x+dx - cx) + square(y - cy);

	if (point1 <= square(r) && point2 <= square(r)) {
		mat[i][j].lx_ratio = 0.0;
	} else if (point1 > square(r) && point2 > square(r)) {
		mat[i][j].lx_ratio = 1.0;
	} else if (point1 <= square(r) && point2 > square(r)) {
		xpoint = cx + sqrt(square(r) - square(y - cy));
		mat[i][j].lx_ratio = (xpoint - x) / dx;
	} else if (point1 > square(r) && point2 <= square(r)) {
		xpoint = cx - sqrt(square(r) - square(y - cy));
		mat[i][j].lx_ratio = (x+dx - xpoint) / dx;
	} else {
		fprintf(stderr, "[Error] Invalid cross section: X-Axis.\n");
		exit(1);
	}
}

// X,Y軸と円柱との交点を解析的に算出(Y軸)
void set_cylinder_ly_ratio(int i, int j, double x, double y, double dx, double dy, double cx, double cy, double r)
{
	double point1, point2;
	double ypoint;

	point1 = square(x - cx) + square(y - cy);
	point2 = square(x - cx) + square(y+dy - cy);

	if (point1 <= square(r) && point2 <= square(r)) {
		mat[i][j].ly_ratio = 0.0;
	} else if (point1 > square(r) && point2 > square(r)) {
		mat[i][j].ly_ratio = 1.0;
	} else if (point1 <= square(r) && point2 > square(r)) {
		ypoint = cy + sqrt(square(r) - square(x - cx));
		mat[i][j].ly_ratio = (ypoint - y) / dy;
	} else if (point1 > square(r) && point2 <= square(r)) {
		ypoint = cy - sqrt(square(r) - square(x - cx));
		mat[i][j].ly_ratio = (y+dy - ypoint) / dy;
	} else {
		fprintf(stderr, "[Error] Invalid cross section: Y-Axis.\n");
		exit(1);
	}
}

void dump_material_matrix()
{
	int i, j;
	FILE *fp1, *fp2, *fp3, *fp4;
	double lx, ly;

	if (bck == 1) {
		fp1 = fopen("wave/material_bck.csv", "w");
	} else {
		fp1 = fopen("wave/material_normal.csv", "w");
	}
	fp2 = fopen("wave/material_area.csv", "w");
	fp3 = fopen("wave/material_lx.csv", "w");
	fp4 = fopen("wave/material_ly.csv", "w");

	lx = 0.0;
	for (i = 0; i < NX; i++) {
		ly = 0.0;
		for (j = 0; j < NY; j++) {
			fprintf(fp1, "%f, %f, %d\n", lx, ly, mat[i][j].master_material);
			fprintf(fp2, "%f, %f, %f\n", lx, ly, mat[i][j].a_ratio);
			fprintf(fp3, "%f, %f, %f\n", lx, ly, mat[i][j].lx_ratio);
			fprintf(fp4, "%f, %f, %f\n", lx, ly, mat[i][j].ly_ratio);
			ly += dy[j];
		}
		fprintf(fp1, "\n"); fprintf(fp2, "\n"); fprintf(fp3, "\n"); fprintf(fp4, "\n");
		lx += dx[i];
	}
	fprintf(fp1, "\n"); fprintf(fp2, "\n"); fprintf(fp3, "\n"); fprintf(fp4, "\n"); 

	fclose(fp1); fclose(fp2); fclose(fp3); fclose(fp4);
}

// FDTDメインルーチン
void calc_fdtd(enum eALG alg)
{
	int i, j;
	int step;
	int last_step = (int)(STOP_TIME/dt);
	double stimulus;
	FILE *fp1, *fp2, *fp3, *fp4, *fp5;

	// Hz[50][50], Hz[90][90]の点のデータを出力するファイル
	if (bck == 1) {
		fp1 = fopen("wave/fdtd_bck.csv", "w");
	} else {
		fp1 = fopen("wave/fdtd_normal.csv", "w");
	}
	print_point_value_header(fp1);
	fp2 = fopen("wave/fdtd_x.csv", "w");
	fp3 = fopen("wave/fdtd_y.csv", "w");
	fp4 = fopen("wave/fdtd_xy.csv", "w");
	fp5 = fopen("wave/fdtd_map.csv", "w");

	// PMLの最外壁を完全電気壁と仮定
	for( i=0; i < total_nx; i++ ) {
		ex[i][0] = 0.0;
		ex[i][total_ny] = 0.0;
	}

	for( j=0; j < total_ny; j++ ) {
		ey[0][j] = 0.0;
		ey[total_nx][j] = 0.0;
	}

	// メインループ
	for (step = 0; step <= last_step; step++) {
		if (step % 100 == 0) {
			printf("%d / %d (%.2f %%)\n", step, last_step, 
						 (double)step/last_step*100);
		}

		// 磁流源の設定
		stimulus = calc_stimulus(step, dt);
		// 波源を(0.01mm, 0.05mm)に設定
		//hz[(int)(0.01/DX)+L][(int)(0.05/DY)+L] += stimulus * dt;
		// 線波源を(x=0.01mm)に設定
		for (j = L; j < NY+L; j++) {
			hz[(int)(0.01/DX)+L][j] += stimulus * dt;
		}

		// Ex&Eyの計算
		std_calc_ex();
		std_calc_ey();

		// PML媒質内の電界計算
		pml_boundary_ex();
		pml_boundary_ey();

		// Hzの計算
		std_calc_hz();

		// PML媒質内の磁界計算
		pml_boundary_hz();

		// 波形出力
		print_point_value(fp1, hz, step*dt, stimulus);
		if (step % 10 == 0) print_line_value(fp2, hz, step, 0);
		if (step % 10 == 0) print_line_value(fp3, hz, step, 1);
//		print_line_value(fp4, hz, step, 2, nx, ny);
		if (step % 10 == 0) print_map_value(fp5, hz, step);
	}

	fclose(fp1); fclose(fp2); fclose(fp3); fclose(fp4);
}

// 入力源
double calc_stimulus(int step, double dt)
{
	//	double TT = 9.4e-9;
	//	double TT = 10e-9;
	double time = step * dt;

	//return square(sin(2.0 * PI * time / TT)) / 2.0;
	return sin(2.0 * PI * time / TT) / 2.0;
}

// 出力ファイルヘッダ
void print_point_value_header(FILE* fp)
{
	fprintf(fp, "#Time[1], Stimulus[2], ");
	fprintf(fp, "P1(0.01,0.05)[3], P2(0.04,0.05)[4], P3(0.06,0.05)[5], P4(0.09,0.05)[6], ");
	fprintf(fp, "P5(0.04,0.01)[7], P6(0.04,0.02)[8], P7(0.04,0.03)[9], P8(0.04,0.04)[10], ");
	fprintf(fp, "P9(0.04,0.06)[11], P10(0.04,0.07)[12], P11(0.04,0.08)[13], P12(0.04,0.09)[14]\n");
}

// データ出力1
void print_point_value(FILE* fp, double** data, double time, double stimulus)
{
	double area = DX * DY;

	fprintf(fp, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n",
					time, stimulus, 
					K * data[(int)(0.01/DX)+L][(int)(0.05/DY)+L] / area, K * data[(int)(0.04/DX)+L][(int)(0.05/DY)+L] / area, 
					K * data[(int)(0.06/DX)+L][(int)(0.05/DY)+L] / area, K * data[(int)(0.09/DX)+L][(int)(0.05/DY)+L] / area,
					K * data[(int)(0.04/DX)+L][(int)(0.01/DY)+L] / area, K * data[(int)(0.04/DX)+L][(int)(0.02/DY)+L] / area,
					K * data[(int)(0.04/DX)+L][(int)(0.03/DY)+L] / area, K * data[(int)(0.04/DX)+L][(int)(0.04/DY)+L] / area,
					K * data[(int)(0.04/DX)+L][(int)(0.06/DY)+L] / area, K * data[(int)(0.04/DX)+L][(int)(0.07/DY)+L] / area,
					K * data[(int)(0.04/DX)+L][(int)(0.08/DY)+L] / area, K * data[(int)(0.04/DX)+L][(int)(0.09/DY)+L] / area);
}

// データ出力2
void print_line_value(FILE* fp, double** data, int step, int mode)
{
	fprintf(fp, "# STEP=%d\n", step);
	int i, min;

	switch(mode) {	
	case 0:   // X軸モード
		for (i = 0; i < total_nx; i++) {
			fprintf(fp, "%d, %g\n", i, data[i][(NY/2)+L]);
//			fprintf(fp, "%d, %g\n", i, data[i][70]);
		}
		break;
	case 1:   // Y軸モード
		for (i = 0; i < total_ny; i++) {
//			fprintf(fp, "%d, %g\n", i, data[90][i]);
			fprintf(fp, "%d, %g\n", i, data[(NX/2)+L][i]);
		}
		break;
	case 2:   // XY軸(45度)モード
		min = (total_nx <= total_ny) ?  total_nx : total_ny;
		for (i = 0; i < min; i++) {
			fprintf(fp, "%d, %g\n", i, data[i][i]);
		}
		break;
	default:
		printf("Error: Unknown axis mode.\n");
		exit(20);
	}
	fprintf(fp, "\n");
}

// データ出力3
void print_map_value(FILE* fp, double** data, int step)
{
	int i, j;
	double lx, ly;

	fprintf(fp, "# STEP = %d\n", step);
	lx = 0.0;
	for (i = 0; i < total_nx; i++) {
		ly = 0.0;
		for (j = 0; j < total_ny; j++) {
			fprintf(fp, "%f, %f, %g\n", lx, ly, data[i][j]);
			ly += dy[j];
		}
		fprintf(fp, "\n");
		lx += dx[i];
	}
	fprintf(fp, "\n");
}

// 配列の初期値を0に設定(DOUBLE)
void initialize_2d_array_dbl(double** ary, int nx, int ny)
{
	int i, j;

	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			ary[i][j] = 0.0;
		}
	}
}

// 配列の初期値を0に設定(CONFORM_T)
void initialize_2d_array_conform_t(conform_t** ary, int nx, int ny)
{
	int i, j;
	
	// 解析領域全体を真空に
	for (i = 0; i < nx; i++) {
		for (j = 0; j < ny; j++) {
			ary[i][j].master_material = 0;
			ary[i][j].second_material = 0;
			ary[i][j].a_ratio = 1.0;
			ary[i][j].lx_ratio = 1.0;
			ary[i][j].ly_ratio = 1.0;
		}
	}
}

// 引数チェック
void check_opt(int argc, char** argv, float* cfl, enum eALG* alg)
{
	if (argc != 5) {
		printf("Error: Invalid Arguments.\n");
		printf("Usage: fdtd cfl NORMAL|ADE|ADE2 BCK_FLAG CALIB_FLAG\n");
		exit(10);
	} else {
		(*cfl) = atof(argv[1]);
		(*alg) = atoi(argv[2]);
		bck = atoi(argv[3]);
		bck_calib = atoi(argv[4]);
	}
}

char* get_bck_mode(int bck)
{
	char* mode;

	if (bck == 1) {
		mode = "ON";
	} else {
		mode = "OFF";
	}

	return mode;
}

char* get_algorithm_name(enum eALG alg)
{
	char* name;

	switch(alg) {
	case NORMAL:
		name = "Normal";
		break;
	case ADE:
		name = "ADE";
		break;
	case ADE2:
		name = "ADE2";
		break;
	default:
		printf("Error: Unkown Method\n");
		exit(12);
	}

	return name;
}

int main(int argc, char** argv)
{
	double min_dx, min_dy;

	float cfl;
	enum eALG alg = ADE;

	// オプションチェック
	check_opt(argc, argv, &cfl, &alg);

	// PML吸収境界層と含めた配列サイズ
	total_nx = NX + 2*L;
	total_ny = NY + 2*L;

	// グリッド配列の確保
	dx = (double *)malloc(sizeof(double) * total_nx);
	dy = (double *)malloc(sizeof(double) * total_ny);

	// グリッド初期化
	set_delta_x();
	set_delta_y();

	// 解析領域の材料
	mat = alloc_2d_array_conform_t(NX, NY);
	initialize_2d_array_conform_t(mat, NX, NY);
	//set_cylinder_material();
	set_rectangle_material();
	dump_material_matrix();

	// CFL条件の設定(BCK込み)
	set_cfl_constant(cfl, &min_dx, &min_dy);
	if (bck == 1) { 
		if (bck_calib == 1) {
			printf("[INFO] Calibrationg Area Ratio of BCK Cells\n");
			calibrate_bck_cell_area(cfl); 
		}
	}

	// 解析情報の出力
	printf("Input Freq : %g [Hz]\n", FREQ);
	printf("Wave Length: %f [m]\n", LIGHT/FREQ);
	printf("Area X : %f [m]\n", PX);
	printf("Area Y : %f [m]\n", PY);
	printf("Grid X : %d\n", NX);
	printf("Grid Y : %d\n", NY);
	printf("Min X  : %f [m]\n", min_dx);
	printf("Min Y  : %f [m]\n", min_dy);
	printf("CFL    : %f\n", cfl);
	printf("dt     : %g [sec]\n", dt);
	printf("Mehtod : %s\n", get_algorithm_name(alg));
	printf("BCK    : %s\n", get_bck_mode(bck));

	// 電磁界配列の確保
	printf("Allocating Memory...");
	alloc_memory_for_analysis_area();

	// 係数初期化
	initialize_coefficient();

	// X方向・Y方向の磁界配列を１次元配列で確保(for PML)
	hzx = (double *)calloc( (total_nx * total_ny) - (NX * NY), sizeof(double) );
	hzy = (double *)calloc( (total_nx * total_ny) - (NX * NY), sizeof(double) );

	printf("Done\n");

	// FDTDメインルーチン
	printf("FDTD Calculating...\n");
	calc_fdtd(alg);
	printf("Finished FDTD Calculating.\n");

	return 0;
}

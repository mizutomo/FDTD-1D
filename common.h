#ifndef __COMMON_H__
#define __COMMON_H__

#define PI            (3.14159)
#define LIGHT         (2.99792458e8)
#define PERMITTIVITY  (8.85418782e-12)
#define PERMEABILITY  (4.0 * PI * 1.0e-7)

// 解析エリアの物理サイズ
#define PX            (1.0)     // 1m

// 解析領域の配列サイズ
#define NX            (1000)

// 解析領域のグリッドサイズ
#define DX            (PX/NX)

// 入射波の周波数
#define FREQ          (100e+6)
#define TT            (1/FREQ)

// CFL条件
#define CFL           (0.5)

// 過渡解析ストップタイム
#define STOP_TIME     (100e-9)

#endif // __COMMON_H__

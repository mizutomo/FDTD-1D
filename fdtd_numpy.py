#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import ctypes
import numpy as np

class Const():
  """ 定数静的クラス """ 
  PI = 3.14159
  LIGHT = 2.99792458e8
  PERMITTIVITY = 8.85418782e-12
  PERMEABILITY = 4.0 * PI * 1.0e-7

  # 解析エリアの物理サイズ
  PX = 1.0     # 1m

  # 解析領域の配列サイズ
  NX = 1000

  # 解析領域のグリッドサイズ
  DX = PX/NX

  # 入射波の周波数
  FREQ = 100e+6
  TT = 1/FREQ

  # CFL条件
  CFL = 0.5

  # 過渡解析ストップタイム
  STOP_TIME = 100e-9
  #STOP_TIME = 10e-9

def square(x):
  """ 二乗お助け関数 """
  return x * x

def find_min_cell(ary):
  """ 最小のセルを見つける関数 """
  return min(ary)

def calc_cfl_constant():
  """ CFL値の算出 """
  return Const.CFL * 1.0 / (Const.LIGHT * (1.0 / Const.DX))

def calc_ey(ey, hz, dx, dt):
  """ CFL値の算出 """
  ey[1:Const.NX] = ey[1:Const.NX] \
      - dt / (Const.PERMITTIVITY * dx[1:Const.NX]) * (hz[1:Const.NX] - hz[0:Const.NX-1])

def calc_hz(hz, ey, dx, dt):
  """ CFL値の算出 """
  hz[0:Const.NX] = hz[0:Const.NX] \
      - dt / (Const.PERMEABILITY * dx[0:Const.NX]) * (ey[1:Const.NX+1] - ey[0:Const.NX])

def ey_boundary(ey, hist_ey, dt):
  """ 境界条件(Mur1st) """
  c = Const.LIGHT * dt
  ey[0] = hist_ey[1] + (c - Const.DX) / (c + Const.DX) * (ey[1] - hist_ey[0])
  ey[Const.NX] = hist_ey[2] + (c - Const.DX) / (c + Const.DX) * (ey[Const.NX-1] - hist_ey[3])

def hist_update(ey, hist_ey):
  """ Mur1st(アップデート) """
  hist_ey[0] = ey[0]
  hist_ey[1] = ey[1]
  hist_ey[2] = ey[Const.NX-1]
  hist_ey[3] = ey[Const.NX]

def print_point_value_header(fp):
  """ データ出力1(ファイルヘッダ) """
  print >> fp, "#Time(1), Stimulus(2), hz[10](3), hz[20](4), hz[30](5), hz[40](6), hz[50](7), hz[60](8), hz[70](9), hz[80](10), hz[90](11)"

def print_point_value(fp, data, time, stimulus):
  """ データ出力1(本体) """
  print >> fp, "%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g" % \
      (time, stimulus, \
         data[100], data[200], data[300], \
         data[400], data[500], data[600], \
         data[700], data[800], data[900])

def print_line_value(fp, data, step):
  """ データ出力2 """
  print >> fp, "# STEP = %d" % step
  for i in range(Const.NX):
    print >> fp, "%d, %g" % (i, data[i])
  print >> fp, ""

# def calc_stimulus(step, dt):
#   """ 入力源(Sin波) """
#   time = step * dt
#   return math.sin(2.0 * Const.PI * time / Const.TT) / 2.0

def calc_stimulus(step, dt):
  """ 入力源(ガウシアンパルス) """
  tc = 10e-9
  pw = 0.1e-9
  t_current = (float(step) - 0.5) * dt
  t_eff = (t_current - tc) / pw;
  return 1.0 * math.exp(-1 * (t_eff * t_eff));

def calc_fdtd(ey, hz, dx, dt, hist_ey):
  """ FDTDメインルーチン"""
  last_step = int(Const.STOP_TIME / dt)

  fp1 = open('fdtd_point.csv', 'w')
  print_point_value_header(fp1)
  fp2 = open('fdtd_line.csv', 'w')

  resource = ctypes.CDLL('utility.dylib')

  tstart = ctypes.c_double()
  resource.get_current_time_by_sec(ctypes.byref(tstart))
  
  # メインループ
  step = 0
  while step <= last_step:
    if step % 100 == 0:
      print "%d / %d (%.2f %%)" % (step, last_step, float(step)/last_step * 100)

    # 磁流源の設定
    stimulus = calc_stimulus(step, dt)
    # 波源を(0.01mm, 0.05mm)に設定
    hz[500] += stimulus * dt

    # Eyの計算
    calc_ey(ey, hz, dx, dt)

    # 境界条件
    ey_boundary(ey, hist_ey, dt)
    hist_update(ey, hist_ey)

    # Hzの計算
    calc_hz(hz, ey, dx, dt)

    # 波形出力
    # print_point_value(fp1, hz, step * dt, stimulus)
    # if step % 100 == 0:
    #   print_line_value(fp2, hz, step)

    step += 1

  tend = ctypes.c_double()
  resource.get_current_time_by_sec(ctypes.byref(tend))
  memsize = ctypes.c_long()
  resource.get_use_memory_size_from_mac(ctypes.byref(memsize))
      
  fp1.close()
  fp2.close()

  print "All User Time: %.2f [sec], %.2f [min], %.2f [hour]" % \
      (tend.value-tstart.value, (tend.value-tstart.value)/60, (tend.value-tstart.value)/3600);

  if 0 <= memsize.value and memsize.value < 1024:
    print "Memory       : %.2f [B]" % (float(memsize.value))
  elif 1024 <= memsize.value and memsize.value/1024 < 1024:
    print "Memory       : %.2f [KB]" % (float(memsize.value)/1024)
  elif 1024 <= memsize.value/1024 and memsize.value/1024/1024 < 1024:
    print "Memory       : %.2f [MB]" % (float(memsize.value)/1024/1024)
  elif 1024 <= memsize.value/1024/1024:
    print "Memory       : %.2f [GB]" % (float(memsize.value)/1024/1024/1024)


def main():
  """ メインルーチン """
  hist_ey = np.array([0.0] * 4)
  
  # グリッド配列の確保&初期化
  dx = np.array([Const.DX] * Const.NX)

  # CFL条件の設定
  dt = calc_cfl_constant()

  # 解析情報の出力
  print "Input Freq : %g [Hz]" % (Const.FREQ);
  print "Wave Length: %f [m]" % (Const.LIGHT/Const.FREQ);
  print "Area X : %f [m]" %  (Const.PX);
  print "Grid X : %d" % (Const.NX);
  print "Min X  : %f [m]" % (Const.DX);
  print "CFL    : %f" % (Const.CFL);
  print "dt     : %g [sec]" % (dt);

  # 電磁界配列の確保
  print "Allocating Memory...", 
  ey = np.array([0.0] * (Const.NX+1))
  hz = np.array([0.0] * Const.NX)
  print "Done"

  # FDTDメインルーチン
  print "FDTD Calculating...";
  calc_fdtd(ey, hz, dx, dt, hist_ey);
  print "Finished FDTD Calculating.";

if __name__ == '__main__':
  main()


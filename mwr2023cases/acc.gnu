#!/usr/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.2 patchlevel 8    last modified 2019-12-01 
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2019
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal qt 0 font "Sans,9"
# set output
unset clip points
set clip one
unset clip two
set errorbars front 1.000000 
set border 31 front lt black linewidth 1.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02 
set style ellipse size graph 0.05, 0.03 angle 0 units xy
set dummy x, y
set format x "% h" 
set format y "% h" 
set format x2 "% h" 
set format y2 "% h" 
set format z "% h" 
set format cb "% h" 
set format r "% h" 
set ttics format "% h"
set timefmt "%d/%m/%y,%H:%M"
set angles radians
set tics back
unset grid
unset raxis
set theta counterclockwise right
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key fixed right top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border  lt -1 linewidth  1.0
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
unset micro
unset minussign
set view 60, 30, 1, 1
set view azimuth 0
set rgbmax 255
set samples 100, 100
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels 5
set cntrparam levels auto
set cntrparam firstlinetype 0 unsorted
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set xyplane relative 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set nomttics
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit autofreq 
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit autofreq 
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit autofreq 
set rtics axis in scale 1,0.5 nomirror norotate  autojustify
set rtics  norangelimit autofreq 
unset ttics
set title "" 
set title  font "" textcolor lt -1 norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" textcolor lt -1 norotate
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ 0.00000 : 200.000 ] noreverse writeback
set x2range [ -279.150 : 331.201 ] noreverse writeback
set ylabel "" 
set ylabel  font "" textcolor lt -1 rotate
set y2label "" 
set y2label  font "" textcolor lt -1 rotate
set yrange [ 0.00000 : 2.50000 ] noreverse writeback
set y2range [ 0.555425 : 2.65595 ] noreverse writeback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse writeback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate
set cbrange [ * : * ] noreverse writeback
set rlabel "" 
set rlabel  font "" textcolor lt -1 norotate
set rrange [ * : * ] noreverse writeback
unset logscale
unset jitter
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "nb_NO.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set pm3d nolighting
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2 size screen 0.05, 0.6 front  noinvert bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
GNUTERM = "qt"
set size 1.0, 1.0
pack( r, g, b ) = 2**16*r + 2**8*g + b
set style line  1  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(179,226,205)  # 41 light green
set style line  2  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(27,158,119)   # 25 dark green
set style line  22 lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(27,158,119)   # 25 dark green
set style line  3  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(244,202,220)  # 44 light red
set style line  4  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(231,41,138)   # 28 dark red
set style line  44 lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(231,41,138)   # 28 dark red
set style line  5  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(203,213,232)   # 43 light blue
set style line  6  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(117,112,179)   # 27 dark blue
set style line  7  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(204,204,204)   # 48 light gray
set style line  8  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(102,102,102)   # 32 dark gray
set style line  9  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(253,205,172)   # 42 light orange
set style line 10  lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(217,95,2)     # 26 dark orange
set autoscale
set yrange [ 0 : 1.5 ] noreverse nowriteback
set xrange [ 0 : 16 ] noreverse nowriteback
set style data linespoints
set nogrid
## Last datafile plotted: "rmse.dat"
set terminal pdfcairo enhanced font  "Arial,40" size 10in,7in lw 2.0 rounded
set key font ",40"

set xlabel "Window lengths"
set ylabel "RMS value"

set output 'rmse-KS-ES.pdf'
p 'rmse-KS-ES.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-ES.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-ES.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-ES.dat' u 1:3 t "rmse Ocean"  linestyle 6

set output 'rmse-KS-ESX.pdf'
p 'rmse-KS-ESX.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-ESX.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-ESX.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-ESX.dat' u 1:3 t "rmse Ocean"  linestyle 6

set output 'rmse-KS-MDA.pdf'
p 'rmse-KS-MDA.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-MDA.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-MDA.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-MDA.dat' u 1:3 t "rmse Ocean"  linestyle 6

set output 'rmse-KS-MDAX.pdf'
p 'rmse-KS-MDAX.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-MDAX.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-MDAX.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-MDAX.dat' u 1:3 t "rmse Ocean"  linestyle 6

set output 'rmse-KS-IES.pdf'
p 'rmse-KS-IES.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-IES.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-IES.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-IES.dat' u 1:3 t "rmse Ocean"  linestyle 6

set output 'rmse-KS-IES4X.pdf'
p 'rmse-KS-IES4X.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-IES4X.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-IES4X.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-IES4X.dat' u 1:3 t "rmse Ocean"  linestyle 6

set output 'rmse-KS-IESX.pdf'
p 'rmse-KS-IESX.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-IESX.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-IESX.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-IESX.dat' u 1:3 t "rmse Ocean"  linestyle 6

set xrange [ 0 : 10 ] noreverse nowriteback
set yrange [ 0 : 1 ] noreverse nowriteback
set xlabel "MDA steps"

set output 'rmse-KS-MDAstep5.pdf'
p 'rmse-KS-MDAstep5.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-MDAstep5.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-MDAstep5.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-MDAstep5.dat' u 1:3 t "rmse Ocean"  linestyle 6



set output 'rmse-KS-MDAstep10.pdf'
p 'rmse-KS-MDAstep10.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-MDAstep10.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-MDAstep10.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-MDAstep10.dat' u 1:3 t "rmse Ocean"  linestyle 6




set output 'rmse-KS-MDAstep15.pdf'
p 'rmse-KS-MDAstep15.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-MDAstep15.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-MDAstep15.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-MDAstep15.dat' u 1:3 t "rmse Ocean"  linestyle 6




set xtics 1.0
set xrange [ 2 : 6 ] noreverse nowriteback
set yrange [ 0 : 0.3 ] noreverse nowriteback
set xlabel "Window lengths"

set output 'rmse-KS-MDAC.pdf'
p 'rmse-KS-MDA.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-MDA.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-MDA.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-MDA.dat' u 1:3 t "rmse Ocean"  linestyle 6

set output 'rmse-KS-IESXC.pdf'
p 'rmse-KS-IESX.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-IESX.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-IESX.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-IESX.dat' u 1:3 t "rmse Ocean"  linestyle 6

set output 'rmse-KS-IES4XC.pdf'
p 'rmse-KS-IES4X.dat' u 1:4 t "rmss Atmos"  linestyle 3,\
  'rmse-KS-IES4X.dat' u 1:5 t "rmss Ocean"  linestyle 5,\
  'rmse-KS-IES4X.dat' u 1:2 t "rmse Atmos"  linestyle 4,\
  'rmse-KS-IES4X.dat' u 1:3 t "rmse Ocean"  linestyle 6

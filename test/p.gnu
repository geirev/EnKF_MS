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
set xrange [ 0.00000 : 1200.00 ] noreverse writeback
set x2range [ 0.00000 : 1023.00 ] noreverse writeback
set ylabel ""
set ylabel  font "" textcolor lt -1 rotate
set y2label ""
set y2label  font "" textcolor lt -1 rotate
set yrange [ -3.00000 : 2.00000 ] noreverse writeback
set y2range [ -2.96401 : 1.83709 ] noreverse writeback
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
x = 0.0
pack( r, g, b ) = 2**16*r + 2**8*g + b
set key noautotitle
set size 1.0, 1.0
set style line  1 lt 1 lw 1 pt 3 ps 0 linecolor rgb pack(179,226,205)  # 41 light green
set style line  2 lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(27,158,119)   # 25 dark green
set style line  22 lt 1 lw 1 pt 3 ps 0 linecolor rgb pack(27,158,119)   # 25 dark green
set style line  3 lt 1 lw 1 pt 3 ps 0 linecolor rgb pack(244,202,220)  # 44 light red
set style line  4 lt 1 lw 4 pt 3 ps 0 linecolor rgb pack(231,41,138)   # 28 dark red
set style line  44 lt 1 lw 1 pt 3 ps 0 linecolor rgb pack(231,41,138)   # 28 dark red
set style line  5 lt 1 lw 3 pt 3 ps 0 linecolor rgb pack(203,213,232)   # 43 light blue
set style line  6 lt 1 lw 3 pt 3 ps 0 linecolor rgb pack(117,112,179)   # 27 dark blue
set style line  7 lt 1 lw 3 pt 3 ps 0 linecolor rgb pack(204,204,204)   # 48 light gray
set style line  8 lt 1 lw 3 pt 3 ps 0 linecolor rgb pack(102,102,102)   # 32 dark gray
set style line  9 lt 1 lw 3 pt 3 ps 0 linecolor rgb pack(253,205,172)   # 42 light orange
set style line 10 lt 1 lw 3 pt 3 ps 0 linecolor rgb pack(217,95,2)     # 26 dark orange
set autoscale
set yrange [ -6 : 6 ] noreverse nowriteback
set xrange [ 0.0 : 1023 ] noreverse nowriteback
set style data linespoints
set nogrid

set xlabel  "distance (x)"

set terminal wxt size 1600,1000

pdfout=0   # Make pdf files
shift=3    # Shifting the y positions of the solutions
cshift=0   # Shifting the y-positions of the covariances
covs=9     # Scaling up the covariance functions
nstd=2     # Number of std dev for uncertainties
f(value, left, right) = (value < left || value > right ? 1/0 : value)

if (pdfout == 1) {set terminal pdfcairo enhanced font  "Arial,15" size 10in,7in lw 1.0 rounded}

set arrow from 341,-6 to 341,6, graph 1 nohead
set arrow from 682,-6 to 682,6, graph 1 nohead

set style fill transparent solid 0.4 noborder
# Initial condition
set title "Time (t=0)"
if (pdfout == 1) {set output "sol_0000I.pdf"}
p 'sol_0000I.dat' u 1:(-shift+$4+2*$6):(-shift+$4-nstd*$6) with filledcurve fc rgb pack(179,226,205) ,\
  'sol_0000I.dat' u 1:(shift+$5+2*$7):(shift+$5-nstd*$7) with filledcurve fc rgb pack(244,202,220) ,\
  'sol_0000I.dat' u 1:(-shift+$2) linestyle 22 title "Ref ocean",\
  'sol_0000I.dat' u 1:(shift+$3) linestyle 44 title "Ref atmos",\
  'sol_0000I.dat' u 1:(-shift+$4) linestyle 2 title "Est ocean",\
  'sol_0000I.dat' u 1:(shift+$5) linestyle 4 title "Est atmos",\
  'sol_0000I.dat' u (f($1,241,441)):(-cshift+covs*$8) linestyle 2 title  "Cov o-o",\
  'sol_0000I.dat' u (f($1,241,441)):(cshift+covs*$9) linestyle 22 title  "Cov  o-a",\
  'sol_0000I.dat' u (f($1,582,782)):(-cshift+covs*$10) linestyle 44 title  "Cov a-o",\
  'sol_0000I.dat' u (f($1,582,782)):(cshift+covs*$11) linestyle 4 title  "Cov a-a"
pause 2.0

do for [var=10:1000:10] {
   xx="000"
   if (var > 9) {xx="00"}
   if (var > 99) {xx="0"}
   if (var > 999) {xx=""}

   list = "F A"
   do for [yy in list] {
      set title "Time (t=".var."".yy.")"
      if (pdfout == 1) {set output "sol_".xx."".var."".yy.".pdf"}
      plot  "sol_".xx."".var."".yy.".dat" u 1:(-shift+$4+nstd*$6):(-shift+$4-nstd*$6) with filledcurve fc rgb pack(179,226,205) ,\
            "sol_".xx."".var."".yy.".dat" u 1:(shift+$5+nstd*$7):(shift+$5-nstd*$7) with filledcurve fc rgb pack(244,202,220) ,\
            "sol_".xx."".var."".yy.".dat" u 1:(-shift+$2) linestyle 22 title  "Ref ocean",\
            "sol_".xx."".var."".yy.".dat" u 1:(shift+$3) linestyle 44 title  "Ref atmos",\
            "sol_".xx."".var."".yy.".dat" u 1:(-shift+$4) linestyle 2 title  "Est ocean",\
            "sol_".xx."".var."".yy.".dat" u 1:(shift+$5) linestyle 4 title  "Est atmos",\
            "oceanobs_".xx."".var.".dat"  u 1:(-shift+$2):(nstd*$3)  w yerr pt 7 ps 1.0 lc rgb pack(27,158,119) title "Obs ocean",\
            "atmosobs_".xx."".var.".dat"  u 1:(shift+$2):(nstd*$3)  w yerr pt 7 ps 1.0 lc rgb pack(231,41,138) title "Obs atmos",\
            "sol_".xx."".var."".yy.".dat" u (f($1,241,441)):(-cshift+covs*$8) linestyle 2 title  "Cov o-o",\
            "sol_".xx."".var."".yy.".dat" u (f($1,241,441)):(cshift+covs*$9) linestyle 22 title  "Cov  o-a",\
            "sol_".xx."".var."".yy.".dat" u (f($1,582,782)):(-cshift+covs*$10) linestyle 44 title  "Cov a-o",\
            "sol_".xx."".var."".yy.".dat" u (f($1,582,782)):(cshift+covs*$11) linestyle 4 title  "Cov a-a"
      pause 0.1
   }
}

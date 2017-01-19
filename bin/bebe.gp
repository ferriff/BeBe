reset

###set boxwidth binwidth
##bw = 5
###bin(x, width) = width * floor(x / width)
##bin(x, width) = width * floor(x / width) + 0.5 * width
##
###plot 'datafile' using (bin($1,bw)):(1.0) smooth freq with boxes
##
##bw = 2
##p 'signal_summary.dat' u (bin($4,bw)):(1.0) sm freq with boxes

pdir = 'plots/'
system("mkdir -p " . pdir)
idir = 'gp_histos'

set terminal pdfcairo enhanced color solid font "Helvetica,18" size 5,4.5

set mxtics 5
set mytics 5
set key samplen 1.5
set tics scale .75
set tmargin 1.5

set auto fix
set offsets graph 0.05, graph 0.05, graph 0.1, graph 0.1
set grid lc rgb '#888888'

save '.setup'


set offsets graph 0.0, graph 0.0, graph 0.1, graph 0.
set out pdir . '/amplitude.pdf'
set xlabel 'Amplitude (ADC counts)'
set ylabel 'Events / 1 ADC count'
set label 1 "x10^{3}" at graph .975, -0.025
p [0:150] idir . '/h_ampl_rebin.dat' u (0.001 * $1):2 not w histeps lt 6 lw 2
set out


set out pdir . '/amplitude_rebin.pdf'
set xlabel 'Amplitude (ADC counts)'
set ylabel 'Events / 16 ADC count'
set label 1 "x10^{3}" at graph .975, -0.025
p [100:106.95] idir . '/h_ampl_rebin.dat' u (0.001 * $1):2 not w histeps lt 6 lw 1.25
set out


set out pdir . '/amplitude_rebin_fit.pdf'
set xlabel 'Amplitude (ADC counts)'
set ylabel 'Events / 16 ADC count'
gauss(x) = N * exp(-(1000.*x-mean)*(1000.*x-mean)/2/sigma/sigma)
#N = 6.07240e+01; mean = 1.01297e+05; sigma = 5.83687e+01;
#N = 5.89721e+01; mean = 1.01279e+05; sigma = 6.11442e+01;
N = 44.2914; mean = 101280; sigma = 60.7363;
set samples 1000
set label 1 "x10^{3}" at graph .975, -0.025
set label 2 sprintf("N = %.1f\n{/Symbol m} = %.0f ADC count\n{/Symbol s} = %.1f ADC count\nFWHM = %.2g%%", N, mean, sigma, 100. * 2.35 * sigma / mean) at graph 0.05, 0.7
p [100.3:101.65] idir . '/h_ampl_rebin.dat' u (0.001 * $1):2 not w histeps lt 6 lw 1.5, gauss(x) not w l lt 7
unset label 1
unset label 2
set out


set out pdir . '/heater_time.pdf'
set xlabel 'Time (s)'
set ylabel 'Events / 0.5 ms'
p [:1500] idir . '/h_dt.dat' u ($1 / 2000.):2 not w histep lw 1.5 lt 6
set out


set out pdir . '/heater_amplitude.pdf'
set xlabel 'Amplitude (ADC counts)'
set ylabel 'Events / 16 ADC count'
set label 1 "x10^{3}" at graph .975, -0.025
p [46.3:46.6] idir . '/h_ampl_heater.dat' u (0.001 * $1):2 not w histep lw 1.5 lt 6, \
              idir . '/h_ampl_heater_drift_corr.dat' u (0.001 * $1):2 not w histep lw 1.5 lt 7
set out


load '.setup'

set out pdir . '/drift.pdf'
set xlabel 'Time (h)'
set ylabel 'Amplitude (ADC count)'
f(x) = m * x + q
m = 1.5; q = 46500
fit f(x) idir . '/g_ampl_vs_t.dat' u 1:2 via m, q
p idir . '/g_ampl_vs_t.dat' u 1:2 not w p pt 7 lt 2 ps .5, f(x) not lt 7 lw 1.5
set out


set out pdir . '/drift2.pdf'
set xlabel 'Time (h)'
set ylabel 'Amplitude (ADC count)'
m = 1.5; q = 101200
fit f(x) idir . '/g_ampl2_vs_t.dat' u 1:2 via m, q
p [][] idir . '/g_ampl2_vs_t.dat' u 1:2 not w p pt 7 lt 6 ps .5, f(x) not lt 7 lw 1.5
set out


set out pdir . '/amplitude_vs_time.pdf'
set xlabel 'Time (h)'
set ylabel 'Amplitude (ADC count)'
p [][0:120000] 'gp_histos/g_ampl_all_vs_t.dat' u 1:2 not w d lt 6
set out


set out pdir . '/amplitude_vs_time_zoom.pdf'
set xlabel 'Time (h)'
set ylabel 'Amplitude (ADC count)'
p [][0:5000] 'gp_histos/g_ampl_all_vs_t.dat' u 1:2 not w d lt 6
set out


set terminal pdfcairo enhanced color solid font "Helvetica,18" size 5,4
set out pdir . '/pedestal.pdf'
set offsets graph 0.0, graph 0.0, graph 0.1, graph 0.0
set ylabel 'Events / ADC count'
set xlabel 'Amplitude (ADC count)'
p [0:3000][] idir . '/h_ped.dat' u 1:2 not w histep lt 6
set out


set out pdir . '/pedestal_rms.pdf'
stat idir . '/h_ped_rms.dat' u 2
ssum = STATS_sum
stat [:15] idir . '/h_ped_rms.dat' u 2
noise = 4.513
set arrow 1 from first noise, graph 0.005 to first noise, graph 0.15
set arrow 1 backhead size 2.5,25 back lw 2 lt 6
set label 10 "{/=14 noise level: 4.513 ADC count}" at first 6, graph 0.1
unset grid
set offsets graph 0.0, graph 0.0, graph 0.1, graph 0.1
set y2tics
set y2range [0:1.05]
set ytics nomirror
set mytics 10
set y2tics auto 0.1
set log y
set y2label 'fraction of Events'
set ylabel 'Events / ADC count'
set xlabel 'Amplitude (ADC count)'
p [0:50] idir . '/h_ped_rms.dat' u 1:($2 / ssum) axis x1y2 not sm cum w l lt 7 dt '.', '' u 1:2 not w histep lt 6
unset log y
set out

reset
l '.setup'

##set out pdir . '/rise_decay_vs_amplitude.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Rise time + decay time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [:150][0:1000] 'summary.dat' u ($1 / 1000):((-$6 + $7) * 2.) not w d lt 6
##set out
##
##
##set out pdir . '/rise_vs_amplitude.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Rise time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [:150][0:150] 'summary.dat' u ($1 / 1000):(($6) * 2.) not w d lt 4
##set out
##
##
##set out pdir . '/decay_vs_amplitude.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Decay time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [:150][0:1000] 'summary.dat' u ($1 / 1000):(($7) * 2.) not w d lt 7
##set out
##
##
##set out pdir . '/rise_decay_inter_vs_amplitude.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Rise time + decay time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [:150][0:1000] 'summary.dat' u ($1 / 1000):((-$8 + $9) * 2.) not w d lt 6
##set out
##
##
##set out pdir . '/rise_inter_vs_amplitude.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Rise time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [:150][0:150] 'summary.dat' u ($1 / 1000):(($8) * 2.) not w d lt 4
##set out
##
##
##set out pdir . '/decay_inter_vs_amplitude.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Decay time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [:150][0:1000] 'summary.dat' u ($1 / 1000):(($9) * 2.) not w d lt 7
##set out
##
##
##set out pdir . '/rise_decay_inter_vs_amplitude_fit.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Rise time + decay time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [0:150][0:800] 'summary.dat' u ($10 / 1000):((-($8 - ($3 - $11)) + $9 + ($3 - $11)) * 2.) not w d lt 6, \
##        '' u ($13 > 10. ? $10 / 1000 : 1/0):((-($8 - ($3 - $11)) + $9 + ($3 - $11)) * 2.) not w d lt 7, \
##        1/0 t 'pedestal r.m.s. <= 10 ADC count' w p pt 7 lt 6, \
##        1/0 t 'pedestal r.m.s.  > 10 ADC count' w p pt 7 lt 7
##set out
##
##
##set out pdir . '/rise_inter_vs_amplitude_fit.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Rise time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [0:150][0:150] 'summary.dat' u ($10 / 1000):(($8 - ($3 - $11)) * 2.) not w d lt 4
##set out
##
##
##set out pdir . '/decay_inter_vs_amplitude_fit.pdf'
##set xlabel 'Amplitude (ADC count)'
##set ylabel 'Decay time (ms)'
##set label 1 "x10^{3}" at graph .975, -0.025
##p [0:150][0:1000] 'summary.dat' u ($10 / 1000):(($9 + ($3 - $11)) * 2.) not w d lt 7
##set out

set macros
set offsets graph 0.0, graph 0.0, graph 0.0, graph 0.

### for 0.2 0.8 thresholds
##triplet = "$13 < 10 && $8 < 19 && $8 > 15 && $10 / 1000. > 3 && $10 / 1000. < 7"
##top1    = "$13 < 10 && $8 > 24 && $8 < 25 && $10 / 1000 < 16"
##top2    = "$13 < 10 && $8 > 23.25 && $8 < 24 && $10 / 1000 > 4"
##top3    = "$13 < 10 && $8 > 22.25 && $8 < 23 && $10 / 1000 < 16"
##top4    = "$13 < 10 && $8 > 21.25 && $8 < 22 + 0.75 * $10 / 100000. && $10 / 1000 > 4"
##top5    = "$13 < 10 && $8 > 19.75 && $8 < 21 && $10 / 1000 > 60"
##heat    = "$13 < 10 && $8 > 19.5  && $8 < 20.5 && $10 / 1000 > 40 && $10 / 1000. < 50"

# for 0.05-0.95 (0.95-0.20) rise (decay) thresholds
triplet = "$13 < 10 && $8 > 30 && $8 < 34 && $10 / 1000. > 3 && $10 / 1000. < 10"
top1    = "$13 < 10 && $8 > 43 && $8 < 44 && $10 / 1000 < 30"
top2    = "$13 < -10 && $8 > 23.25 && $8 < 24 && $10 / 1000 > 4"
top3    = "$13 < 10 && $8 > 41.1 + 0.02 * $10 / 1000. && $8 < 42.5 + 0.01 * $10 / 1000. && $10 / 1000 < 60"
top4    = "$13 < 10 && $8 > 40.25 + 0.0075 * $10 / 1000.  && $8 < 40.85 + 0.015 * $10 / 1000. && $10 / 1000 > 15"
top5    = "$13 < 10 && $8 > 38.75 + 0.01 * $10 / 1000. && $8 < 39.5 + .01 * $10 / 1000. && $10 / 1000 > 10"
top6    = "$13 < 10 && $8 > 38 && $8 < 39 && $10 / 1000 > 90"
heat    = "$13 < 10 && $8 > 36  && $8 < 39 && $10 / 1000 > 40 && $10 / 1000. < 50"

#bkg  = 
#sig  = 
heat = "$13 < 10 && $8 > 36  && $8 < 39 && $10 / 1000 > 40 && $10 / 1000. < 50"

set out pdir . '/colour_rise_vs_amplitude.pdf'
set xlabel 'Amplitude (ADC count)'
set ylabel 'Rise time (ms)'
set label 1 "x10^{3}" at graph .975, -0.025
p [0:140][18.5:21.5] 'summary.dat' u ($13 < 10 ? $10 / 1000 : 1/0):($8 * 0.5) not w d lt 6, \
         '' u (@heat ? $10 / 1000. : 1/0):($8 * 0.5) not w d lt 7
#p [0:140][18.5:21.5] 'summary.dat' u ($13 < 10 ? $10 / 1000 : 1/0):($8 * 0.5) not w d lc rgb '#cccccc', \
#         '' u (@top1 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top2 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top3 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top4 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top5 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@heat ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@triplet ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top6 ? $10 / 1000. : 1/0):($8 * 0.5) not w d
set out


set out pdir . '/colour_rise_vs_amplitude_zoom.pdf'
set xlabel 'Amplitude (ADC count)'
set ylabel 'Rise time (ms)'
set label 1 "x10^{3}" at graph .975, -0.025
p [0:140][19.5:21] 'summary.dat' u ($13 < 10 ? $10 / 1000 : 1/0):($8 * 0.5) not w d lt 6
#p [0:140][19.5:21] 'summary.dat' u ($13 < 10 ? $10 / 1000 : 1/0):($8 * 0.5) not w d lc rgb '#cccccc', \
#         '' u (@top1 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top2 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top3 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top4 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top5 ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@heat ? $10 / 1000. : 1/0):($8 * 0.5) not w d bt 6, \
#         '' u (@triplet ? $10 / 1000. : 1/0):($8 * 0.5) not w d, \
#         '' u (@top6 ? $10 / 1000. : 1/0):($8 * 0.5) not w d
set out


set out pdir . '/colour_decay_vs_amplitude.pdf'
set xlabel 'Amplitude (ADC count)'
set ylabel 'Decay time (ms)'
set label 1 "x10^{3}" at graph .975, -0.025
p [0:140][90:130] 'summary.dat' u ($10 / 1000):($9 * 0.5) not w d lt 6, \
         '' u (@heat ? $10 / 1000. : 1/0):($9 * 0.5) not w d lt 7
#p [0:140][90:210] 'summary.dat' u ($10 / 1000):($9 * 0.5) not w d lc rgb '#cccccc', \
#         '' u (@top1 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top2 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top3 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top4 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top5 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@heat ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@triplet ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top6 ? $10 / 1000. : 1/0):($9 * 0.5) not w d
set out


set out pdir . '/colour_decay_vs_amplitude_zoom.pdf'
set xlabel 'Amplitude (ADC count)'
set ylabel 'Decay time (ms)'
set label 1 "x10^{3}" at graph .975, -0.025
p [0:140][95:112] 'summary.dat' u ($10 / 1000):($9 * 0.5) not w d lt 6, \
         '' u (@heat ? $10 / 1000. : 1/0):($9 * 0.5) not w d lt 7
#p [0:140][100:125] 'summary.dat' u ($10 / 1000):($9 * 0.5) not w d lc rgb '#cccccc', \
#         '' u (@top1 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top2 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top3 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top4 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top5 ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@heat ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@triplet ? $10 / 1000. : 1/0):($9 * 0.5) not w d, \
#         '' u (@top6 ? $10 / 1000. : 1/0):($9 * 0.5) not w d
unset label 1
set out

exit

#set out pdir . '/pulses.pdf'
#set ylabel 'a.u.'
#set xlabel 'Time (ms)'
#flags = system("grep '^#' pulses.dat | tr -d '#flag: ' | tr '\n' ' '")
#print sprintf("%d pulse(s) found.", words(flags))
#p [290:400][0:1.1] for[i=2:words(flags):1] 'pulses.dat' u (0. + word(flags, i) > 0 ? ($1 + $2) * 0.5 : 1/0):($3) index i - 1 not w d lt word(flags, i)
#set out


set out pdir . '/pulses_zoom_rise.pdf'
set ylabel 'a.u.'
set xlabel 'Time (ms)'
flags = system("grep '^#' pulses.dat | tr -d '#flag: ' | sort | uniq | tr '\n' ' '")
do for [i=1:words(flags)] {
        flag = word(flags, i)
        print "selecting pulses with flag ".flag."..."
        system("grep '^#flag: ".flag."' -A 2002 pulses.dat | sed -e '/--/d' > pulses_flag_".flag.".dat")
}
#print sprintf("%d pulse(s) found.", words(flags))
#p [190:230][0:1.1] for[i=2:words(flags):1] 'pulses.dat' u (0. + word(flags, i) > 0 ? ($1 + $2) * 0.5 : 1/0):($3) index i - 1 not w d lt word(flags, i)
#set out


set offsets graph 0.0, graph 0.0, graph 0.1, graph 0.1

set out pdir . '/one_pulse.pdf'
set xlabel 'sample number'
set ylabel 'a.u.'
set ytics auto .25
set arrow 10 from first 180, .95 to first 275, .95 nohead lt 7
set arrow 11 from first 170, .05 to first 250, .05 nohead lt 7
set arrow 12 from first 300, .20 to first 400, .20 nohead lt 7
#set label 1 "x10^{3}" at graph -0.05, 1.1
stats 'gp_one_pulse/Graph.dat' u ($2 * 0.001)
M = STATS_max
#p [-400:][-.4:] 'gp_one_pulse/Graph.dat' u ($1 / 2.):($2 * 0.001 / M) not w d lt 6
p [][] 'gp_one_pulse/Graph.dat' u ($1 / 2.):($2 * 0.001 / M) not w d lt 6
set out


set out pdir . '/one_pulse_maxfit.pdf'
set xlabel 'sample number'
set ylabel 'a.u.'
#set label 1 "x10^{3}" at graph -0.05, 1.1
stats 'gp_one_pulse/Graph.dat' u ($2 * 0.001)
M = STATS_max
farabola(x) = p0 + p1 * x + p2 * x * x
parabola(x) = x > 453 && x < 463 ? p0 + p1 * x + p2 * x * x : 1/0
fit [453:463] farabola(x) 'gp_one_pulse/Graph.dat' u ($1):($2 * 0.001 / M) via p0,p1,p2
p [445:473] 'gp_one_pulse/Graph.dat' u ($1):($2 * 0.001 / M) not w p pt 7 ps .25 lt 6, parabola(x) not lt 7
set out
#
#
#set out pdir . '/one_signal_rise_decay.pdf'
#set xlabel 'sample number'
#set ylabel 'a.u.'
##set label 1 "x10^{3}" at graph -0.05, 1.1
#stats 'gp_one_pulse/Graph.dat' u ($2 * 0.001)
#M = STATS_max
#p [400:1000][0.1:0.3] 'gp_one_pulse/Graph.dat' u ($1):($2 * 0.001 / M) not w p pt 7 ps .25 lt 6, parabola(x) not lt 7

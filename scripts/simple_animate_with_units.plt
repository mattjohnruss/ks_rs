load "parula_lines.pal"

if(exist("gif_output")) {
    set term gif size 1920,1080 animate delay 3
    set output "animate.gif"
    delay = 0.0
}
else {
    set term qt noraise font "Sans,12"
}

set clip two

set key top right

#fn(i) = sprintf("< sed '/^$/d' output_%05i.csv", i)
fn(i) = sprintf("output_%05i.csv", i)
fn_exact(i) = sprintf("output_exact_%05i.csv", i)

n_var = words(vars)

# Avogadro's number
n_a = 6.02214076e23

# number of minutes per unit time
time_scale = 2500.0 / 60.0

# number of metres per unit space
space_scale = 5.0e-4

# number of (mol/litre) per unit concentration for chemokines
chemokine_scale = 1.0e-8
#chemokine_scale = 1.0e-8 * n_a

# number of (mol/litre) per unit concentration for cells
cell_scale = 6.40825486e-14
#cell_scale = 6.40825486e-14 * n_a

# numbers of (mol/litre) per unit concentration for each variable
array var_scales[7]
var_scales[1] = chemokine_scale
var_scales[2] = chemokine_scale
var_scales[3] = chemokine_scale
var_scales[4] = cell_scale
var_scales[5] = cell_scale
var_scales[6] = cell_scale
var_scales[7] = cell_scale

set xlabel "x (m)"
set ylabel "Concentration (mol/litre)"
#set ylabel "Concentration (molecules or cells/litre)"

if(exists("n_inc")) {
}
else {
    n_inc = 1
}

if(exists("n_start")) {
}
else {
    n_start = 0
}

do for [i = n_start:n:n_inc] {
    stats [*:*] [*:*] fn(i) u 1 skip 1 nooutput
    set label 1 sprintf("t = %.1f min (%.3f non-dim.)", time_scale*STATS_min, STATS_min) at graph 0.05,0.95

    if(exist("steady")) {
        plot for [j = 1:n_var] fn(i)               u 2:2+word(vars,j) w l lc j lw 2 ti columnheader, \
             for [j = 1:n_var] "output_steady.csv" u 2:2+word(vars,j) w l lc j lw 2 dt 3 noti
    }
    else {
        if(exist("exact")) {
            plot for [j = 1:n_var] fn(i)           u 2:2+word(vars,j) w l lc j lw 2 ti columnheader, \
                 for [j = 1:n_var] fn(i)           u 2:(exactf($2))   w l lc j dt 2 ti "exact"
        }
        else {
            if(exist("exact_file")) {
                plot for [j = 1:n_var] fn(i)       u 2:2+word(vars,j) w l lc j lw 1 ti columnheader, \
                     for [j = 1:n_var] fn_exact(i) u 2:2+word(vars,j) w l lc j lw 2 dt 3 ti columnheader
            }
            else {
                plot for [j = 1:n_var] fn(i)       u (space_scale*$2):(var_scales[word(vars,j) + 0]*column(2+word(vars,j))) w l lc j lw 2 ti col(2+word(vars,j))
            }
        }
    }

    if(exist("delay")) {
        pause(delay)
    }
    else {
        pause("0.02")
    }
}

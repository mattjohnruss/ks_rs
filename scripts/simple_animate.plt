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

fn(i) = sprintf("< sed '/^$/d' output_%05i.csv", i)
#fn(i) = sprintf("output_%05i.csv", i)
fn_exact(i) = sprintf("output_exact_%05i.csv", i)

n_var = words(vars)

set xlabel "x"
set ylabel "Concentration"

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
    set label 1 sprintf("t = %.3f", STATS_min) at graph 0.05,0.95

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
                plot for [j = 1:n_var] fn(i)           u 2:2+word(vars,j) w l lc j lw 1 ti columnheader, \
                     for [j = 1:n_var] fn_exact(i)     u 2:2+word(vars,j) w l lc j lw 2 dt 3 ti columnheader
            }
            else {
                plot for [j = 1:n_var] fn(i)           u 2:2+word(vars,j) w l lc j lw 2 smooth unique ti columnheader
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

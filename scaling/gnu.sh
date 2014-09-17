set logscale
set terminal png
set out "strong.png"
set title "Strong Scaling (16-2048 cores)"
set xlabel "# cores"
set ylabel "Speedup (16 cores = 16)"
set key bottom right
set xrange [16:2048]
set yrange [16:2048]
plot "scaling_4_8.txt" u 1:(68.3*16/$2) w linespoints t "4 levels (4681 sub-grids)", "scaling_5_8.txt" u 1:(619/$2*16) w linespoints t "5 levels (37449 sub-grids)", x t "Perfect Scaling"

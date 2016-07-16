set key bottom right

set xlabel "distance to channel center [nm]"
set ylabel "all quantities in simulation units"

plot \
  "eof_electrokinetics.dat" u 1:2  w p pt 6 lc rgb "#ff0000" ti "EK density",\
  "eof_analytical.dat" u 1:2  w l      lc rgb "#ff8080" ti "analytic density",\
  "eof_electrokinetics.dat" u 1:3  w p pt 6 lc rgb "#0000ff" ti "LB velocity",\
  "eof_analytical.dat" u 1:3  w l      lc rgb "#8080ff" ti "analytic velocity",\
  "eof_electrokinetics.dat" u 1:4  w p pt 6 lc rgb "#00a000" ti "LB stress xz",\
  "eof_analytical.dat" u 1:4  w l      lc rgb "#50a050" ti "analytic stress xz",\

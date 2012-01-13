set term post eps enh "Times-Roman" 25
set out "nacl-rdf.eps"
set xlabel "r"
set ylabel "g(r)"
plot [.9:0.5*6.75] \
	"data/rdf_from_melt_00.data" notitle w linesp pt 4, \
	"data/rdf_from_melt_10.data" notitle w linesp pt 6, \
	"data/rdf_lj_00.data"        notitle w linesp pt 1
unset out

set out "neutral-rho.eps"
set xlabel "z"
set ylabel "{/Symbol r}(z)"
plot [0.8:3.5] \
	"data/neutral-rho.data" u 1:2 notitle w linesp pt 4, \
	"data/neutral-rho.data" u 1:3 notitle w linesp pt 5
unset out

set out "nonneutral-rho.eps"
set xlabel "z"
set ylabel "{/Symbol r}(z)"
plot [0.8:3.5] \
	"data/nonneutral-rho.data" u 1:2 notitle w linesp pt 4, \
	"data/nonneutral-rho.data" u 1:3 notitle w linesp pt 5
unset out

set cnt 0
for {set i 0} {$i < 100} {incr i} { lappend avg_rdf 0 }

foreach filename [lrange $argv 1 end] {

    set f [open $filename "r"]
    while { [blockfile $f read auto] != "eof" } {}
    close $f

    set rdf [analyze rdf 0 1 0.9 [expr $box_l/2] 100]
    set rlist ""
    set rdflist ""
    foreach value [lindex $rdf 1] {
	lappend rlist [lindex $value 0]
	lappend rdflist [lindex $value 1]
    }
    
    set avg_rdf [vecadd $avg_rdf $rdflist]
    incr cnt
}

set avg_rdf [vecscale [expr 1.0/$cnt] $avg_rdf]

set plot [open "rdf.data" "w"]
puts $plot "\# r rdf(r)"
foreach r $rlist rdf $avg_rdf { puts $plot "$r $rdf" }
close $plot

set n_part 200; set density 0.7
set box_l [expr pow($n_part/$density,1./3.)]

setmd box_l $box_l $box_l $box_l
setmd periodic 1 1 1

set q 1; set type 0
for {set i 0} { $i < $n_part } {incr i} {
  set posx [expr $box_l*[t_random]]
  set posy [expr $box_l*[t_random]]
  set posz [expr $box_l*[t_random]]
  set q [expr -$q]; set type [expr 1-$type]
  part $i pos $posx $posy $posz q $q type $type
}

setmd time_step 0.01; setmd skin 0.3
set temp 1; set gamma 1
thermostat langevin $temp $gamma

set sig 1.0; set cut [expr 1.12246*$sig]
set eps 1.0; set shift [expr 0.25*$eps]
inter 0 0 lennard-jones $eps $sig $cut $shift 0
inter 1 0 lennard-jones $eps $sig $cut $shift 0
inter 1 1 lennard-jones $eps $sig $cut $shift 0
puts [inter coulomb 10.0 p3m tunev2 accuracy 1e-3 mesh 32]

set p3m_params [inter coulomb]
foreach f $p3m_params { eval inter $f }

if { [regexp "ROTATION" [code_info]] } {
    set deg_free 6
} {
    set deg_free 3
}

set integ_steps 200
for {set cap 20} {$cap < 200} {incr cap 20} {
  puts "t=[setmd time] E=[analyze energy total]"
  inter ljforcecap $cap
  integrate $integ_steps
}
inter ljforcecap 0

for {set i 0} { $i < 20 } { incr i} {
    set temp [expr [analyze energy kinetic]/(($deg_free/2.0)*$n_part)]
    puts "t=[setmd time] E=[analyze energy total], T=$temp"
    integrate $integ_steps

    set f [open "config_$i" "w"]
    blockfile $f write tclvariable {box_l density}
    blockfile $f write variable box_l
    blockfile $f write particles {id pos type}
    close $f
}

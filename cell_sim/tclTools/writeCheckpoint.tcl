proc writeCheckpoint { filename } {

  set f [open "$filename.cpt" "w"]
  blockfile $f write tclvariable all
  blockfile $f write variable all
  blockfile $f write particles {id pos v virtual} all
  blockfile $f write interactions
  blockfile $f write bonds all
  close $f
  
  lbfluid save_ascii_checkpoint "$filename.cptlb"
}
import ctypes
import sys
sys.setdlopenflags((sys.getdlopenflags() | ctypes.RTLD_GLOBAL ))

import espresso as es
import numpy

# this should be replaced by a python call in future:
es._espressoHandle.Tcl_Eval("thermostat langevin 1. 1.")

N=10
VarId=0;

print "\nTest global variables: set them (if not RO) and compare the Tcl and python output:\n"

es.glob.skin=1.
varname="skin";
# set the global variable 
es.glob.skin=0.01
# get it
py_val=es.glob.skin;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd skin"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="box_l";
# set the global variable 
es.glob.box_l=[10., 10., 10.]
# get it
py_val=es.glob.box_l;
# check if it is accessible to the original Tcl interface
tcl_val=es._espressoHandle.Tcl_Eval("setmd box_l").split();
tcl_val=numpy.array([ float(tcl_val[0]), float(tcl_val[1]), float(tcl_val[2]) ]);
VarId=VarId+1;
for i in range(py_val.size):
  if py_val[i] != tcl_val[i]:
    raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="cell_grid";
# get it
py_val=es.glob.cell_grid;
# check if it is accessible to the original Tcl interface
tcl_val=es._espressoHandle.Tcl_Eval("setmd cell_grid").split();
tcl_val=numpy.array([ float(tcl_val[0]), float(tcl_val[1]), float(tcl_val[2]) ]);
VarId=VarId+1;
for i in range(py_val.size):
  if py_val[i] != tcl_val[i]:
    raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="cell_size";
# get it
py_val=es.glob.cell_size;
# check if it is accessible to the original Tcl interface
tcl_val=es._espressoHandle.Tcl_Eval("setmd cell_size").split();
tcl_val=numpy.array([ float(tcl_val[0]), float(tcl_val[1]), float(tcl_val[2]) ]);
VarId=VarId+1;
for i in range(py_val.size):
  if py_val[i] != tcl_val[i]:
    raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="dpd_gamma";
# get it
py_val=es.glob.dpd_gamma;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd dpd_gamma"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="dpd_r_cut";
# get it
py_val=es.glob.dpd_r_cut;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd dpd_r_cut"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="gamma";
# get it
py_val=es.glob.gamma;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd gamma"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="integ_switch";
# get it
py_val=es.glob.integ_switch;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd integ_switch"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="local_box_l";
# get it
py_val=es.glob.local_box_l;
# check if it is accessible to the original Tcl interface
tcl_val=es._espressoHandle.Tcl_Eval("setmd local_box_l").split();
VarId=VarId+1;
tcl_val=numpy.array([ float(tcl_val[0]), float(tcl_val[1]), float(tcl_val[2]) ]);
for i in range(py_val.size):
  if py_val[i] != tcl_val[i]:
    raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="max_cut";
# get it
py_val=es.glob.max_cut;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd max_cut"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="max_num_cells";
# set the global variable 
es.glob.max_num_cells=125
# get it
py_val=es.glob.max_num_cells;
# check if it is accessible to the original Tcl interface
tcl_val=int(es._espressoHandle.Tcl_Eval("setmd max_num_cells"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="min_num_cells";
# set the global variable 
es.glob.min_num_cells=16
# get it
py_val=es.glob.min_num_cells;
# check if it is accessible to the original Tcl interface
tcl_val=int(es._espressoHandle.Tcl_Eval("setmd min_num_cells"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="max_part";
# get it
py_val=es.glob.max_part;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd max_part"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="max_range";
# get it
py_val=es.glob.max_range;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd max_range"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="max_skin";
# get it
py_val=es.glob.max_skin;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd max_skin"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="n_layers";
# get it
py_val=es.glob.n_layers;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd n_layers"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="n_nodes";
# get it
py_val=es.glob.n_nodes;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd n_nodes"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="n_part";
# get it
py_val=es.glob.n_part;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd n_part"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="n_part_types";
# get it
py_val=es.glob.n_part_types;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd n_part_types"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="n_rigidbonds";
# get it
py_val=es.glob.n_rigidbonds;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd n_rigidbonds"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="node_grid";
# set the global variable 
es.glob.node_grid=[1., 1., 1.]
# get it
py_val=es.glob.node_grid;
# check if it is accessible to the original Tcl interface
tcl_val=es._espressoHandle.Tcl_Eval("setmd node_grid").split();
tcl_val=numpy.array([ float(tcl_val[0]), float(tcl_val[1]), float(tcl_val[2]) ]);
VarId=VarId+1;
for i in range(py_val.size):
  if py_val[i] != tcl_val[i]:
    raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="nptiso_gamma0";
# get it
py_val=es.glob.nptiso_gamma0;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd nptiso_gamma0"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="nptiso_gammav";
# get it
py_val=es.glob.nptiso_gammav;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd nptiso_gammav"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="npt_p_ext";
# get it
py_val=es.glob.npt_p_ext;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd npt_p_ext"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="npt_p_inst";
# get it
py_val=es.glob.npt_p_inst;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd npt_p_inst"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="npt_p_inst_av";
# get it
py_val=es.glob.npt_p_inst_av;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd npt_p_inst_av"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="npt_p_diff";
# set it
es.glob.npt_p_diff=12.0;
# get it
py_val=es.glob.npt_p_diff;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd npt_p_diff"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="npt_piston";
# set it
es.glob.npt_piston=11.0;
# get it
py_val=es.glob.npt_piston;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd npt_piston"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="periodicity";
# set the global variable 
es.glob.periodicity=[1., 1., 1.]
# get it
py_val=es.glob.periodicity;
# check if it is accessible to the original Tcl interface
tcl_val=es._espressoHandle.Tcl_Eval("setmd periodicity").split();
tcl_val=numpy.array([ float(tcl_val[0]), float(tcl_val[1]), float(tcl_val[2]) ]);
VarId=VarId+1;
for i in range(py_val.size):
  if py_val[i] != tcl_val[i]:
    raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="temperature";
# get it
py_val=es.glob.temperature;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd temperature"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="thermo_switch";
# get it
py_val=es.glob.thermo_switch;
# check if it is accessible to the original Tcl interface
tcl_val=int(es._espressoHandle.Tcl_Eval("setmd thermo_switch"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="time";
# set it
es.glob.time=131.130;
# get it
py_val=es.glob.time;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd time"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="time_step";
# set the global variable 
es.glob.time_step=0.01
# get it
py_val=es.glob.time_step;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd time_step"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";
  
varname="timings";
# set the global variable 
es.glob.timings=17
# get it
py_val=es.glob.timings;
# check if it is accessible to the original Tcl interface
tcl_val=es._espressoHandle.Tcl_Eval("setmd timings");
if str(py_val) != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";
  
varname="transfer_rate";
# get it
py_val=es.glob.transfer_rate;
# check if it is accessible to the original Tcl interface
tcl_val=int(es._espressoHandle.Tcl_Eval("setmd transfer_rate"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="max_cut_nonbonded";
# get it
py_val=es.glob.max_cut_nonbonded;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd max_cut_nonbonded"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


varname="verlet_reuse";
# get it
py_val=es.glob.verlet_reuse;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd verlet_reuse"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="lattice_switch";
# get it
py_val=es.glob.lattice_switch;
# check if it is accessible to the original Tcl interface
tcl_val=int(es._espressoHandle.Tcl_Eval("setmd lattice_switch"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="dpd_tgamma";
# get it
py_val=es.glob.dpd_tgamma;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd dpd_tgamma"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="dpd_tr_cut";
# get it
py_val=es.glob.dpd_tr_cut;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd dpd_tr_cut"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="dpd_twf";
# get it
py_val=es.glob.dpd_twf;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd dpd_twf"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="dpd_wf";
# get it
py_val=es.glob.dpd_wf;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd dpd_wf"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="adress_vars";
# get it
py_val=es.glob.adress_vars;
# check if it is accessible to the original Tcl interface
tcl_val=es._espressoHandle.Tcl_Eval("setmd adress_vars").split();
tcl_val=numpy.array([ \
  float(tcl_val[0]), float(tcl_val[1]), 
  float(tcl_val[2]), \
  float(tcl_val[3]), \
  float(tcl_val[4]), \
  float(tcl_val[5]), \
  float(tcl_val[6]) ]);
VarId=VarId+1;
for i in range(py_val.size):
  if py_val[i] != tcl_val[i]:
    raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";

varname="max_cut_bonded";
# get it
py_val=es.glob.max_cut_bonded;
# check if it is accessible to the original Tcl interface
tcl_val=float(es._espressoHandle.Tcl_Eval("setmd max_cut_bonded"));
VarId=VarId+1;
if py_val != tcl_val:
  raise ValueError(varname + " FAILED\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");
print (str(VarId)+" "+varname).ljust(20), "OK";


# print the last varuable if desired
# print("\n" + varname + ":\n" + "Tcl".ljust(10) + str(tcl_val) + "\n" +  "python".ljust(10) + str(py_val) + "\n");

print "Everything OK"

es._espressoHandle.die()


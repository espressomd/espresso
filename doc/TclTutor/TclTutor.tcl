#!/bin/sh 
# \
exec wish "$0" "$@"

;# NAME:   TclTutor.tcl
;# AUTHOR: Clif Flynt
;# DATE:   Apr 22, 2000
;# DESC:   
;#         
;#
;# PARAMETERS:  Tutor indices can be set from cmd line.
;#
;# RCSID: $Header$
;#
;# This code copyright 2000 to Clif Flynt
;# 9300 Fleming Rd
;# Dexter, MI  48130
;# clif@cflynt.com
;# All Rights Reserved.


# console show

# Clear existing options - we'll set what we want.
option clear

# Delete existing windows - while debugging from tkcon.tcl.

foreach ch [winfo children .] {
    destroy $ch
}

#
# Check of obsolete Wish interpreters
#

if {[catch {namespace eval foo {}}]} {
    text .t  -height 4
    .t insert 0.0 "Your Wish Interpreter ($tk_version) is too old.\n You need to get 8.0 or more recent from dev.scriptics.com"
    grid .t
    button .b -text quit -command exit
    grid .b -row 1 -column 0
    vwait foo
}

######################################################################
# I need these procs to do the setup.  All other procs wait until
#   initial setup is complete to be defined.
######################################################################

################################################################
# proc checkWrap {}--
#    Determines if the application is wrapped, using Dennis Labelle's 
# freeWrap or Jan Nijtmans' wraper.
# Arguments
#   NONE
# 
# Results
#   return "wrap"   if wrapped with Nijtman-wrap
#   return "free"   if wrapped with freeWrap
#   return "" if not wrapped.

proc checkWrap {} {
  global argv0
  if {![string match [info command ::wrap::open] ""]} {
    return "wrap"
  }

  if {(![string match [info command _freewrap_iswrapped] ""])
        && ( ![string match [set f [_freewrap_iswrapped $argv0]] ""])} {
    close $f
    return "free"
  }
return ""
}

################################################################
# proc Debug {arg}--
#    Display debug messages if Tutor(debug) is != 0
# Arguments
#   A string to display
# 
# Results
#   A string is printed to stdout on Unix boxes, and put into a 
#   scrolling textwidget in a toplevel on Windows/Mac platforms
# 
proc Debug {arg} {
    global Tutor tcl_platform

    if {$Tutor(debug) == 0} {return}

    switch $tcl_platform(platform) {
      "unix"	{
          puts "$arg"
      }
      "windows" 	-
      "mac"	{
          if {(![info exists Tutor(debugWin)]) || 
	      (![winfo exists $Tutor(debugWin)])} {
		set r 0
                toplevel .debug -background $Tutor(opt.Frame.background)
		set t .debug.txt

                text $t -yscrollcommand "${t}_scroll set"
                scrollbar ${t}_scroll -command "$t yview"
                grid $t -row $r -column 0 -sticky news 
                grid ${t}_scroll -row $r -column 10 -sticky ns
                grid rowconfigure . $r -weight 10
                grid columnconfigure . $r -weight 10
                set Tutor(debugWin) $t
	  }
      $Tutor(debugWin) insert end "$arg\n"
      }
      default	{
          puts "WHAT: NO $tcl_platform(platform)"
      }
   }
}

################################################################
# proc parseArgs {stateVar {throwError 1}}--
#    Parses $argv into a state array.
#    looks for arguments of the form -key val1 ?val2 ... valn?
#    and assigns them as ArrayName(key) "val1 ?val2 ... valn?"
#     Keys must have a default value to be assigned.
#     By default, an error is thrown if an undefaulted key is found.
#     If throwError == 0, then the undefaulted keays are appended
#        into a list of errors, which will be returned.
# Arguments
#   stateVar    The name of the associative array to have key/values set
#   throwError  Optional 0 to turn off errors and return a list.
# 
# Results
#   Any default values in stateVar are replaced with new values from the
#       command line.
# 

proc parseArgs {stateVar {throwError 1}} {
    global argv
    upvar $stateVar stateArray

    set errs ""
    
    foreach arg $argv {
        if {![string first "-" $arg]} {  
            set index [string range $arg 1 end]
            if {![info exists stateArray($index)]} {
                if {$throwError} {
                    error "No default for ${stateVar}($index)"
                } else {
                    lappend errs "No default for ${stateVar}($index)"
                }
            }
            set cmd set
        } else {
            if {[info exists cmd]} {
                $cmd stateArray($index) $arg
                set cmd lappend
            }
        }
    }
    return $errs
}


# Set platform dependant paths.

    switch $tcl_platform(platform) {
      "unix"	{
            set Tutor(sourceHome) [file dirname $argv0]
            set Tutor(lessonHome) [file dirname $argv0]
	    set Tutor(rcHome) $env(HOME)
            set Tutor(rcfile) [file join $Tutor(rcHome) .tcltutorrc]
            set Tutor(logFileName) [file join $Tutor(rcHome) .tcltutoract]
	    set Tutor(fontMod) -4
       }
      "windows" 	{
            set Tutor(sourceHome) [file attributes [file dirname $argv0] -shortname]
            set Tutor(lessonHome) [file attributes [file dirname $argv0] -shortname]
	    set Tutor(rcHome) [file dirname $argv0]
            set Tutor(rcfile) [file join $Tutor(rcHome) tcltutor.rc]
            set Tutor(logFileName) [file join $Tutor(rcHome) tcltutor.act]
	    set Tutor(fontMod) -5
      }
      "mac"	{
            set Tutor(sourceHome) [file dirname $argv0]
            set Tutor(lessonHome) [file dirname $argv0]
	    set Tutor(rcHome) [file dirname $argv0]
            set Tutor(rcfile) [file join $Tutor(rcHome) tcltutor.rc]
            set Tutor(logFileName) [file join $Tutor(rcHome) tcltutor.act]
	    set Tutor(fontMod) -5
       }
      default	{
          puts "WHAT: NO Support for: $tcl_platform(platform)"
	  exit;
      }
   }

#
# Override the platform specifics if this is wrapped with Jan Nijtmans'
#  wrap, or Dennis Labelle's Freewrap.
#

if {[string match [checkWrap] "wrap"]} {
    set Tutor(wrapNamespace) "::wrap::"
    set Tutor(sourceHome) ""
    set Tutor(lessonHome) ""
    set Tutor(Tcl.lessonFile) Tcl0.lsn
}

if {[string match [checkWrap] "free"]} {
    set Tutor(sourceHome) "G:/TCL_STUFF/TclTutor/TclTutor"
    set Tutor(lessonHome) "G:/TCL_STUFF/TclTutor/TclTutor"
    set Tutor(Tcl.lessonFile) "G:/TCL_STUFF/TclTutor/TclTutor/Tcl0.lsn"
}


# Initialize the various state indices.

set Tutor(version) "2.0 beta4"
set Tutor(debug) 0

set Tutor(courseName) Tcl
set Tutor(courseLevel) 1
set Tutor(Tcl.lessonFile) [file join $Tutor(lessonHome) Tcl0.lsn]
set Tutor(geometry) "500x400"
set Tutor(unique) 0
set Tutor(interpreter) interp
set Tutor(lessonInfo) "Welcome to TclTutor :: Click Help"
set Tutor(fontSize) $Tutor(fontMod)

set Tutor(logFile) ""
set Tutor(logUsage) 0
set Tutor(mailUsage) 0
set Tutor(logCount) 0

set Tutor(terseness.0) Expert
set Tutor(terseness.1) User
set Tutor(terseness.2) Beginner

set Tutor(errorBackground)	     yellow
set Tutor(opt.Button.background)     #38e
set Tutor(opt.Menubutton.background) #38e
set Tutor(opt.Frame.background)      #bff
set Tutor(opt.Label.background)      #bff
set Tutor(opt.Canvas.background)     #bff
set Tutor(opt.Scrollbar.background)  #bff
set Tutor(opt.Text.background)       #9df

set Tutor(opt.Button.font) {helvetica 10}
set Tutor(opt.Menubutton.font) {helvetica 10}
set Tutor(opt.Label.font) {helvetica 10}

set Tutor(opt.Button.foreground)     #000
set Tutor(opt.Menubutton.foreground) #000
set Tutor(opt.Frame.foreground)      #000
set Tutor(opt.Label.foreground)      #000
set Tutor(opt.Canvas.foreground)     #000
set Tutor(opt.Scrollbar.foreground)  #000
set Tutor(opt.Text.foreground)       #000

set Tutor(opt.Button.highlightbackground)     #000
set Tutor(opt.Menubutton.highlightbackground) #000
set Tutor(opt.Frame.highlightbackground)      #000
set Tutor(opt.Label.highlightbackground)      #000
set Tutor(opt.Canvas.highlightbackground)     #bff
set Tutor(opt.Scrollbar.highlightbackground)  #000
set Tutor(opt.Text.highlightbackground)       #000
set Tutor(wrapNamespace) ""

# These group widgets that need the same foreground/background values
#  into single units.

set Tutor(grp.Buttons) {Button Menubutton}
set Tutor(grp.Labels) {Frame Canvas Label Scrollbar}
set Tutor(grp.Text) {Text}

wm title . "TclTutor   $Tutor(version)"




# Parse once to set home, and rcfile if necessary

parseArgs Tutor

catch {$Tutor(wrapNamespace)source $Tutor(rcfile)}

# Redo parse to let command line overrrule settings from rcfile

parseArgs Tutor

# Load the script libraries.
#  Using packages would be good, but would require some more installation
#  overhead.  This will work without making things complex.

foreach f [list scaler.tcl htmllib.tcl] {
    eval $Tutor(wrapNamespace)source [file join $Tutor(sourceHome) $f]
}

# Get the color for frames, and configure the primary window for that color.

frame .f
set Tutor(background) [option get .f background Frame ]
. configure -background $Tutor(background)
destroy .f

#
# Local definitions of procedures win over those in libraries.
#


################################################################
# proc getLessonTitles {}--
#    Searches through the lesson directory for lessons that match the
#    current lesson type, and finds the titles for each lesson.
# Arguments
#   None
# 
# Results
#   Returns a list of lesson number, name and title.

proc getLessonTitles {} {
  global Tutor

  if {[string match [checkWrap] ""]} {
    set filelist [glob \
      [file join $Tutor(lessonHome) $Tutor(courseName)\[0-9\]\[0-9\].lsn] \
      [file join $Tutor(lessonHome) $Tutor(courseName)\[0-9\].lsn] ]
  } else {
    set filelist "G:/TCL_STUFF/TclTutor/TclTutor/Tcl0.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl1.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl10.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl11.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl12.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl13.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl14.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl15.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl16.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl17.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl18.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl19.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl2.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl20.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl21.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl22.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl23.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl24.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl25.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl26.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl27.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl28.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl29.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl3.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl30.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl31.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl32.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl33.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl34.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl35.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl36.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl37.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl38.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl39.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl4.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl40.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl41.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl42.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl43.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl5.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl6.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl7.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl8.lsn G:/TCL_STUFF/TclTutor/TclTutor/Tcl9.lsn"
  }

  foreach file $filelist {
    set infl [$Tutor(wrapNamespace)open $file ]
    set gottitle 0;
    while {!$gottitle} {
      set len [gets $infl line]

      if {$len > 2} {
        if {[string first ":TITLE:" $line] >= 0} {
	  set gottitle 1;
          regsub ":TITLE:" $line "" line
	  set Tutor(lsn.title) $line
          regsub ".*$Tutor(courseName)(\[0-9\]*).lsn" $file {\1} num
	  lappend titlelist [list $num $file $Tutor(lsn.title)]
	  }
	}
      }
    close $infl;
    }

  set Tutor(lessonList) [lsort -command cmpLessonLists $titlelist]

  return $titlelist;
}

################################################################
# proc getCourseTitles {}--
#    Return a list of available courses for creating the menu
# Arguments
#   NONE
# 
# Results
#   No Side effects
# 
proc getCourseTitles {} {
  global Tutor

  if {[string match [checkWrap] ""]} {
    set filelist [glob [file join $Tutor(lessonHome) *.cfg]]
  } else {
    set filelist "Tcl"
  }

  foreach file $filelist {
        set f [file tail $file]
        regsub {.cfg} $f "" courseName
	lappend titlelist [list "" $file $courseName]
    }
  set Tutor(courseList) $titlelist;
  return $titlelist;
} 

################################################################
# proc count {start end {incr 1}}--
#    Return a list of numbers from start to end, 
# Arguments
#   start	The first value to be returned
#   end		The last value will be less than this
#   incr	What to increment by - defaults to 1
# 
# Results
#   No Side Effects
# 
proc count {start end {incr 1}} {
    for {set i $start} {$i < $end} {incr i $incr} {
        lappend rslt $i
    }
    return $rslt
}



################################################################
# proc FillMenu {menu list cmd}--
#    Fills a cascading menu from a list of lessons or courses
# Arguments
#   menu	The menu widget to have items added to it
#   list	A list of courses/lessons
#   cmd		The command to evaluate when a command menu is selected.
# Results
#   

proc FillMenu {menu list cmd} {
  global Tutor
  
  destroy $menu 
  menu $menu

  set result 1;
  
  
  set length [llength $list]
  for {set i 0} {$i < $length} {incr i 10} {
     set last [expr $i + 10];
     if {$last > $length} {set last $length}

     if {$length > 10} {
        $menu add cascade -label "Lesson $i - [expr $last-1]" -menu $menu.lst$i
	set use [menu $menu.lst$i]
	} else {
	  set use $menu
     }

     for {set j $i} {$j < $last} {incr j} {
        set lesson [lindex $list $j]
        foreach {num file title} $lesson {}
	if {[string length $num] > 0} {
	    $use add command -label "$num: $title" -command "$cmd $file"
	} else {
	    $use add command -label "$title" -command "$cmd $file"
	}
     }
  }
}

################################################################
# proc selectCourse {cfgFileName}--
#    The procedure evaluated when a course is selected from the
#    "Select Course" menu
# Arguments
#   cfgFileName		The name of the configuration file to source.
# 
# Results
#   A new config file is sourced, the current lesson is loaded for 
#   that course.
# 
proc selectCourse {cfgFileName} {
    global Tutor
    $Tutor(wrapNamespace)source $cfgFileName
    set Tutor(courseName) [file root [file tail $cfgFileName]]
    getLessonTitles
    FillMenu .file.mnu.lessons $Tutor(lessonList) showLesson

    if {![info exists Tutor($Tutor(courseName).lessonFile)]} {
        set Tutor($Tutor(courseName).lessonFile) \
	    [file join $Tutor(lessonHome) "$Tutor(courseName)0.lsn"]
    }

    showLesson $Tutor($Tutor(courseName).lessonFile)
}

################################################################
# proc cmpLessonLists {a b}--
#    Compares two lesson lists - used to sort the lists
# Arguments
#   a, b	Two lesson lists.
# 
# Results
#   
# 
proc cmpLessonLists {a b} {
    return [expr [lindex $a 0] > [lindex $b 0]]
}

################################################################
# proc readLesson {lessonFilename}--
#    Read and parse a lesson.
# Arguments
#   lessonFileName	The name of the lesson file to read
# 
# Results
#   Tutor indices lsn* (0, 1, 2, code, setup) are initialized from
#   the file. 
# 
proc readLesson {lessonFilename} {
    global Tutor tcl_platform

  set fail [catch {$Tutor(wrapNamespace)open $lessonFilename} infl]
  if {$fail} {
      global errorInfo
      Debug "Can't open $lessonFilename"
      Debug $errorInfo
  }
  foreach t [array names Tutor lsn.*] {
      set Tutor($t) ""
  }
  
  set Tutor(lsn.codeMod) ""
  set Tutor(lsn.setup) ""
  
  while {[set len [gets $infl line]] >= 0} {

      set save 1

      if {[regexp {:LESSON_TEXT_START_LEVEL ([0-9]):} $line m1 level]} {
          regsub {:LESSON_TEXT_START_LEVEL [0-9]:} $line "" line
	  set dest "lsn.$level"
	  set save 0
      }
      if {[regsub {:TEXT_END:} $line "" line]} {
          set dest "none"
	  set save 0
      }

      if {[regsub {:CODE_START:} $line "" line]} {
          set dest "lsn.code"
	  set save 0
	  eval $Tutor(lsn.setup)
      }

      if {[regsub {RCSID:} $line "RCSID:" line]} {
          set dest "none"
	  set save 0
      }

      if {[regsub {::CMD::} $line "" line]} {
          set dest "lsn.setup"
      }

      if {[regsub ":TITLE:" $line "" line]} {
          set dest "none"
      }
      
      if {$save} {
          if {[string match $dest "lsn.code"]} {
              foreach cmd $Tutor(lsn.codeMod) {
	          eval $cmd
	      }
	  }
          append Tutor($dest) "$line\n"
      }
      
  }

  close $infl

}

################################################################
# proc showLesson {lessonFilename}--
#    Loads a lesson file (invoking readLesson) and displays the
#    contents of the appropriate level, code, etc.
# Arguments
#   lessonFilename	The name of the lessonfile to display
# 
# Results
#   The display is updated to reflect the new lesson.
#   Tutor(lsn*) is updated.
#   Tutor(lessonInfo) is updated to reflect new lesson
# 
proc showLesson {lessonFilename} {
    global Tutor

    set Tutor($Tutor(courseName).lessonFile) $lessonFilename

    set list [lindex $Tutor(lessonList) [getLessonIndex]]
    foreach {num file title} $list {}
    set Tutor(lessonInfo) "#$num: $title   ---   Level: $Tutor(terseness.$Tutor(courseLevel))"
    
    logUsage "Lesson: $Tutor(courseName) [file tail $lessonFilename] Verbosity: $Tutor(courseLevel)"
    
    readLesson $lessonFilename
    
    $Tutor(text.lesson) configure -state normal

    foreach b [list runexample nextlesson previouslesson] {
        $Tutor(button.$b) configure -state disabled
    }

    foreach t [list .lesson .code .output] {
        $Tutor(text$t) delete 0.0 end
    }

    HMparse_html $Tutor(lsn.$Tutor(courseLevel)) "HMrender $Tutor(text.lesson) "
    
    $Tutor(text.code) insert 0.0 $Tutor(lsn.code)

    $Tutor(text.lesson) configure -state disabled

    foreach b [list runexample nextlesson previouslesson] {
        $Tutor(button.$b) configure -state normal
    }
    
}

################################################################
# proc setOptions {}--
#    Finds all the State variables associated with widget options, and
#    does the appropriate "option add "  commands.
# Arguments
#   No Options
# 
# Results
#   The default foreground/background colors are set using the option,
#   command, so the widget creating commands don't need tons of 
#   -foo bar settings.
# 
proc setOptions {} {
    global Tutor

    foreach i [array names Tutor opt*] {
        regsub "opt." $i "*" id
        option add $id $Tutor($i)
    }
}


################################################################
# proc createDisplay {}--
#    Generates the 3 window display
# Arguments
#   None
# 
# Results
#   The primary window gets tons of widgets and windows.
# 
proc createDisplay {} {
    global Tutor
    global Scalers
    global tcl_platform
    
    setOptions

    . configure -background $Tutor(opt.Frame.background)

    # Set the rows to low values for scaling - will reset them with scaler
    
    for {set i 0} {$i <9} {incr i} {
        grid rowconfigure . $i -weight 0
        grid columnconfigure . $i -weight 1
    }
    
    # Create the text windows and scrollbars for lessons, code and output
    
    foreach t [list .lesson .code .output] r [count 2 8 2] {
        text $t -yscrollcommand "${t}_scroll set"  -height 1 
        scrollbar ${t}_scroll -command "$t yview"
        grid $t -row $r -column 0 -sticky news -columnspan 9
        grid ${t}_scroll -row $r -column 10 -sticky ns
        grid rowconfigure . $r -weight 10
        set Scalers($r) 10
	set Tutor(text$t) $t
        bind $Tutor(text$t) <Double-Button-1> {+
			showManEntry %W %x %y
			};

    }
    
    HMinit_win .lesson
    
    # Add scalers for the text windows.
    
    grid [makeScaler . 2 $Tutor(opt.Frame.background) "Lesson"] -row 3 -column 0 -columnspan 9
    grid [makeScaler . 4 $Tutor(opt.Frame.background) "Example Code"] -row 5 -column 0 -columnspan 9
    grid [makeScaler . 6 $Tutor(opt.Frame.background) "Output"] -row 7 -column 0 -columnspan 9
    
    menubutton .file -text "File" -menu .file.mnu -border 3 -relief raised
    menu .file.mnu
    grid .file -row 0 -column 1

    
    getLessonTitles

    .file.mnu add cascade -label "Lessons" -menu .file.mnu.lessons 
    FillMenu .file.mnu.lessons $Tutor(lessonList) showLesson
    
    .file.mnu add cascade -label "Courses" -menu .file.mnu.courses
    FillMenu .file.mnu.courses [lsort [getCourseTitles]] selectCourse
    
    foreach l [list "Set Font Size" "Set Colors" "Exit"] {
        regsub -all " " $l "" l2
        .file.mnu add command -label $l -command $l2
    }
    
    if {[string match $tcl_platform(platform) "unix"]} {
        if {$Tutor(logUsage)} {
           .file.mnu add command -label "Disable activity log" -command {set Tutor(logUsage) 0}
         } else {
           .file.mnu add command -label "Enable activity log" -command {set Tutor(logUsage) 1}
	}

        if {$Tutor(mailUsage)} {
           .file.mnu add command -label "Disable mailing log" -command {set Tutor(mailUsage) 0}
	 } else {
           .file.mnu add command -label "Enable mailing log" -command {set Tutor(mailUsage) 1}
	}
    }

    menubutton .terse -text "Terseness" -menu .terse.mnu -border 3 -relief raised
    menu .terse.mnu
    grid .terse -row 0 -column 2
    
    foreach n [count 0 3] t [list Beginner "User" "Expert"] {
        .terse.mnu add command -label $t -command "set Tutor(courseLevel) [expr 2 - $n]; showLesson \$Tutor(\$Tutor(courseName).lessonFile)"
    }


    foreach col [count 3 6] l [list "Next Lesson" "Previous Lesson" "Run Example" ] {
      regsub " " $l  "" c
      set b .[string tolower $c]
      button $b -text $l -command $c -highlightbackground $Tutor(opt.Frame.background)
      grid $b -row 0 -column $col
      set Tutor(button$b) $b
    }

    menubutton .help -text "Help/About" -menu .help.mnu -border 3 -relief raised
    menu .help.mnu
    grid .help -row 0 -column 6
    
    foreach n [count 0 2] t [list Help About] {
        .help.mnu add command -label $t -command "displayHtmlFile \
	    [file join $Tutor(lessonHome) [string tolower $t].html]"
    }

    
    label .lessoninfo -textvar Tutor(lessonInfo) -width 80
    grid .lessoninfo -row 1 -column 0 -columnspan 9
    wm geometry . $Tutor(geometry)
    
    setFonts
}

################################################################
# proc SaveState {}--
#    Saves the contents of the Tutor state array var.
# Arguments
#   NONE
# 
# Results
#   An rcfile is opened, and new data put into it.
# 
proc SaveState {} {
    global Tutor

    # Remember where and how large this window was
    set Tutor(geometry) [winfo geometry .]

    # Clear the log file - this channel won't be here when we start again.
    set Tutor(logFile) ""
    
    # And dump the state information.
    set of [open $Tutor(rcfile) "w"]
    puts $of "array set Tutor [list [array get Tutor]]"
    close $of
}

################################################################
# proc Exit {}--
#    Blow this popsicle stand
# Arguments
#   none
# 
# Results
#   Saves state, and exits
# 
proc Exit {} {
  logUsage "EXIT"
  SaveState
  exit
}

################################################################
# proc getLessonIndex {}--
#    Finds a lesson by file name, and 
#    returns the index of a lesson in the lesson list
# Arguments
#   NONE
# 
# Results
#   No side effects
# 
proc getLessonIndex {} {
    global Tutor

    set index [lsearch $Tutor(lessonList) "*[file tail $Tutor($Tutor(courseName).lessonFile)]*"]
    return $index 
}

################################################################
# proc NextLesson {}--
#    Load the next lesson
# Arguments
#   none
# 
# Results
#   Loads the next lesson.  All Tutor(lsn*) are updated.
# 
proc NextLesson {} {
    global Tutor
    
    set index [getLessonIndex]

    incr index;
    if {$index >= [llength $Tutor(lessonList)]} {
        incr index -1
    }
    set l [lindex $Tutor(lessonList) $index]

    set Tutor($Tutor(courseName).lessonFile) [lindex $l 1]
    showLesson [lindex $l 1]
}

################################################################
# proc PreviousLesson {}--
#    Load the Previous lesson
# Arguments
#   none
# 
# Results
#   Loads the Previous lesson.  All Tutor(lsn*) are updated.
# 
proc PreviousLesson {} {
    global Tutor
    
    set index [getLessonIndex]

    incr index -1;
    if {$index < 0} {
        incr index 1
    }
    set l [lindex $Tutor(lessonList) $index]

    set Tutor($Tutor(courseName).lessonFile) [lindex $l 1]
    showLesson [lindex $l 1]
}




################################################################
# proc dummyputs {args}--
#    A puts for the slave interpreter on Windows or Macs
# Arguments
#   args	The arguments that you'd give to puts
# 
# Results
#   If 'normal' puts, copies the output to global var "output"
#   If going to a pipe, it sends the output to the pipe via 'realputs'

proc dummyputs {args} {
  global Tutor child

#   append Tutor(example.Output) "CMD IS: $args"

  if {([llength $args] > 3) || ([llength $args] < 1)} {
    error "bad argument : should be puts ?-nonewline? ?stream? text \NOT: $args"
    }

   switch "[llength $args]" {
   {1} {
       set args [lindex $args 0]
       append Tutor(example.Output) "$args\n"  
       #   To debug in standalone mode.
       # $child eval [list realputs $args]
       }
   {2} {
       if {[string match "-nonewline" [lindex $args 0]]} {
         set args [lindex $args 1]
         append Tutor(example.Output) "$args"
       } else {
         $child eval realputs $args
        }
    }   
   {3} {
       if {([string match "-nonewline" [lindex $args 0]]) ||
           ([string match "nonewline" [lindex $args 2]])} {
         $child eval realputs $args
        } else {
        error "bad argument : should be puts ?-nonewline? ?stream? text \NOT: $args"
        }
   }
   default {
        error "DFLT: bad argument : should be puts ?-nonewline? ?stream? text \NOT: $args"
   }
 }  
}   

################################################################
# proc RunExample {}--
#    Runs the code in the example code window.
#    Reports and displays errors as best it can.
# Arguments
#   NONE
# 
# Results
#   New data is displayed in Tutor(lsn.output).
#   code window may have a highlighted section if there are errors.
# 

proc RunExample {} {
    global Tutor argv argv0 argc tk_version errorInfo child
    $Tutor(text.output) delete 0.0 end;
    set Tutor(example.Status) 0;
    set Tutor(example.Output) ""

    $Tutor(text.code) tag delete BAD

    if {![string match $Tutor(interpreter) "interp"]} {
       runExternal
       return
    } 

  set txt [subst "set argv0 \"$argv0\";\n"]
  append txt [subst "set argv \"$argv\";\n"]
  append txt [subst "set argc $argc;\n"]

  append txt "global argv argv0 env argc;\n"
  append txt [$Tutor(text.code) get 0.0 end]
  
  set child [interp create]

  if {[string match [string tolower $Tutor(courseName)] "tk"]} {
    load $Tutor(libTk) tk $child
    }
  

  $child eval rename puts realputs
  $child alias puts dummyputs

proc killChildInterp.$child {} " \
      update idle; \
      $child eval update idle; \
      interp delete $child; "

  $child alias exit "killChildInterp.$child"

  set errorInfo ""
  
  set cmd ""

  foreach l [split $txt "\n"] {  
    append cmd "$l\n"
    if {[info complete $cmd]} {
      set fail [evalNdisplay $cmd $child]
      set Tutor(example.Output) ""
      set cmd ""
      if {$fail} {
          break;
	  }
    }
  }

  if {![string match $cmd ""]} {
    set fail [evalNdisplay $cmd $child]
  }

  if {$fail} {
      logUsage "Run Code: error"
  } else {
      logUsage "Run Code: OK"
  }
}

################################################################
# proc evalNdisplay {cmd child}--
#    Evaluate a chunk of code, and display output.
#    If an error occurs, parse the error output
# Arguments
#   cmd		The Tcl command to evaluate.
#   child	The child interpreter
# 
# Results
#   May change the display if there is an error or generated output.
# 
proc evalNdisplay {cmd child} {
    global Tutor argv argv0 argc tk_version errorInfo

      set fail [catch {$child eval $cmd} result]
      ShowOutput "$Tutor(example.Output)"
      if {$fail} {
	  set rlist [ParseError [$Tutor(text.code) get 0.0 end ] $cmd $result $errorInfo]
	  eval $Tutor(text.code) tag add BAD $rlist
          $Tutor(text.code) tag configure BAD -background $Tutor(errorBackground)
	  Debug "ERROR: RLIST $rlist"
	  regsub {invoked from within
"$child eval $cmd"} $errorInfo "" errorInfo
	  ShowOutput "\n--------\n$errorInfo"
	  return 1
	  } else {
	  }
      return 0
}

################################################################
# proc ShowOutput {list}--
#    Display a list of lines
# Arguments
#   list	A list of lines to display in the output window.
# 
# Results
#   New data is displayed.
# 
proc ShowOutput {list} {
  global Tutor
#  set lines [split $list "\n"]
#puts "SHOWOUTP: $list"
#  foreach line $lines {
#         set xx [string trim $line]
#         $Tutor(text.output) insert end "$xx\n"
#    }
         $Tutor(text.output) insert end "$list"

  }

################################################################
# proc setWidgetColor {g opt}--
#    Set a color option for a specific type of widget, button, label, etc
# Arguments
#   g		The index for this option in Tutor (opt.Button.background)
#   opt		The option to set (background, foreground, etc.)
# 
# Results
#   The Tutor value associated with this widget/option is modified.
#   The screen is redrawn to reflect the new value.
# 
proc setWidgetColor {g opt} {
  global Tutor
  
  set wid [lindex [split $g "."] 1]

  set n [format "opt.%s.%s" [lindex $Tutor($g) 0] $opt]
  
  set color [tk_chooseColor -initialcolor $Tutor($n) -title "Set Color for $wid $opt" ]

  if {![string match $color ""]} {
      foreach w $Tutor($g) {
          set n [format "opt.%s.%s" $w $opt]
	  set Tutor($n) $color
	  if {[string match $opt "background"]} {
              set n [format "opt.%s.%s" $w highlightbackground]
	      set Tutor($n) $color
	  }
      }
      configureRecursive .
  }
    setOptions
}

################################################################
# proc SetColors {}--
#    Invoked from the Set Colors menu choice.
#    This creates the toplevel widget for selecting foreground/background colors
# Arguments
#   NONE
# 
# Results
#   May result in new colors being selected.
# 
proc SetColors {} {
  global Tutor
  
  set t [toplevel .colorset_$Tutor(unique)  -background $Tutor(opt.Frame.background)]
  
  incr Tutor(unique)

  set col 0;
  
  foreach g [lsort [array names Tutor grp*]] {
      set wid [lindex [split $g "."] 1]

      label $t.l$col -text $wid
      grid $t.l$col -row 1 -col $col
      set row 2
      foreach opt {foreground background} {
         button $t.b$wid$opt -text $opt -command "setWidgetColor $g $opt"
         grid $t.b$wid$opt -col $col -row $row
	 incr row
      }
      incr col
  }

  set b [button $t.quit -text "Done" -command "setOptions; destroy $t"]
  grid $b -row 0 -column 2

}

################################################################
# proc configureRecursive {parent}--
#    Goes through all the windows and sets new foreground/background
#    values
# Arguments
#   parent	A parent window or frame
# 
# Results
#   Display is redrawn to reflect state variables.
# 
proc configureRecursive {parent} {
  global Tutor

  $parent configure -background $Tutor(opt.Frame.background)
  $parent configure -highlightbackground $Tutor(opt.Frame.highlightbackground)

  foreach ch [winfo children $parent]  {
     switch [set id [winfo class $ch]] {
	 Toplevel -
	 Frame {
	     configureRecursive $ch
	 }
	 Canvas -
         Scrollbar -
	 Button -
	 Text -
	 Label -
	 Scale -
	 Menubutton {
	   foreach n [array names Tutor opt.$id.*] {
	       foreach {o t option} [split $n "."] {}
	       catch {$ch configure -$option $Tutor($n)}
	   }
	 }
	 default {
	   Debug "OOPS, not handling $id"
	 }
     }
  }
}

################################################################
# proc applySize {}--
#    Apply the new font size request.
# Arguments
#   None: Reads value from scaler widget textvariable
# 
# Results
#   Screen is redrawn with new sized fonts.
# 
proc applySize {} {
    global Tutor

    set Tutor(fontSize) [expr $Tutor(tmp) + $Tutor(fontMod)]

    setFonts

    showLesson $Tutor($Tutor(courseName).lessonFile)
}

################################################################
# proc setFonts {}--
#    Set the fonts in the html widget, and text windows based on
#    current value in Tutor(fontSize).
# Arguments
#   NONE
# 
# Results
#   New default fonts.
# 
proc setFonts {} {
    global Tutor HM_globals

    set HM_globals(S_adjust_size) $Tutor(fontSize)
    $Tutor(text.output) configure  -font [HMx_font Courier 14 medium r]
    $Tutor(text.code)   configure  -font [HMx_font Courier 14 medium r]

}

################################################################
# proc SetFontSize {}--
#    Creates a toplevel window for selecting font sizes.
# Arguments
#   NONE
# 
# Results
#   May redraw screen with new sized letters.
# 
proc SetFontSize {} {
  global Tutor
  set t [toplevel .sizeSet_$Tutor(unique) -background $Tutor(opt.Frame.background) ]
  incr Tutor(unique)
  
  wm title $t "Font Size"
  
  set Tutor(tmp) [expr $Tutor(fontSize) - $Tutor(fontMod)]
  
  scale $t.sc -from 0 -to 20 -length 100 -showvalue 1 \
      -variable Tutor(tmp) -orient horizontal \
      -background $Tutor(opt.Frame.background)
  
  
  
  button $t.quit -text "Cancel" -command "destroy $t"
  button $t.apply -text "Apply" -command "applySize"
  button $t.done -text "Done" -command "applySize; destroy $t"
  
  grid $t.sc -row 0 -column 0 -columnspan 3
  grid $t.quit -row 1 -column 0 
  grid $t.apply -row 1 -column 1
  grid $t.done -row 1 -column 2
  
}


################################################################
# proc logUsage {str}--
#    Dump a usage line to the activity log.
#    If more than 10 lines in log, mail it.
# Arguments
#   str		The string to place in the log
# 
# Results
# If not previously opened, the log is opened/created.
# The log is larger.  
# 
proc logUsage {str} {
    global Tutor
    
    if {$Tutor(logUsage) == 0} {return}

    if {[string match $Tutor(logFile) ""]} {
        set Tutor(logFile) [open $Tutor(logFileName) "a"]
    }
   
    set tm [clock format [clock seconds] -format "%d/%b/%y %H:%M:%S"]
    puts $Tutor(logFile) "$tm : $str"
    incr Tutor(logCount)
    if {($Tutor(mailUsage)) && ($Tutor(logCount) > 10)} {
        mailLog
	set Tutor(logCount) 0;
    }
}

################################################################
# proc mailLog {}--
#    Ship it- send the activity log to me
# Arguments
#   None
# 
# Results
#   The log hits the e-mail system, and gets emptied.
# 
proc mailLog {} {
  global Tutor tcl_platform
  catch {close $Tutor(logFile)}
  
    switch $tcl_platform(platform) {
      "unix"	{
	    exec mail tutoractivity@cflynt.com < $Tutor(logFileName)
       }
      "windows" 	-
      "mac"	{

       }
      default	{
          puts "WHAT: NO Support for: $tcl_platform(platform)"
	  exit;
      }
   }
  set Tutor(logFile) [open $Tutor(logFileName) "w"]
}

################################################################
# proc htmlWindow {text}--
#    Open a toplevel window with a scrollbar and an HTML text widget
#    Display the text in that window.
# Arguments
#   text	An HTML document to display.
# 
# Results
#   A new toplevel widget is created.  Unique is incremented.
# 
proc htmlWindow {text} {
    global Tutor
    
    set win .w_$Tutor(unique);
    incr Tutor(unique)
    toplevel $win  -background $Tutor(opt.Frame.background)
    set r 0
    
    set t [text $win.t -yscrollcommand "${win}.scroll set"]
    set sc [scrollbar ${win}.scroll -command "$t yview"]

    grid $t -row $r -column 0 -sticky news
    grid $sc -row $r -column 1 -sticky ns
    grid rowconfigure $win $r -weight 10
    grid columnconfigure $win 0 -weight 10

    HMinit_win $t

    HMparse_html $text "HMrender $t"
    
    button $win.b -text "DONE" -command "destroy $win"
    grid $win.b -row 2 -column 0
}

################################################################
# proc displayHtmlFile {root}--
#    Display a file in a toplevel window.
# Arguments
#   fileName 	The root of the html file name.
# 
# Results
#   
# 
proc displayHtmlFile {fileName} {
    global Tutor
   set infl [$Tutor(wrapNamespace)open $fileName "r"]
   set data [read $infl]
   close $infl
   htmlWindow $data
}

################################################################
# proc showManEntry {window xloc yloc}--
#    Displays a man page for the word at xloc,yloc in the text window.
#    If necessary, opens a new application for this.
# Arguments
#   window	The name of the text widget
#   xloc	X pixel location of a part of the word
#   yloc	Y pixel location of a part of the word
# 
# Results
#   A new window will be created, or an running man page reader
#   will have the display updated to reflect the selected word.

proc showManEntry {window xloc yloc} {
  global Tutor tcl_platform

  set manWord [$window get "@$xloc,$yloc wordstart" "@$xloc,$yloc wordend"]

  if {[info exists tcl_platform(platform)] && [string match "windows" $tcl_platform(platform)]} {
    
    # Get the help file name - several steps of directory manipulations
    #  and taking the last item in the possible list.
    set hlp [file attributes [file dirname [info nameofexecutable]] -shortname]
    set hlp [lindex [glob [file join [file dirname $hlp] doc *.hlp]] end]

    if {[file exists $hlp ]} {
      if {[string first $manWord [info commands]] >= 0} {
        exec C:/Windows/winhlp32.exe -k$manWord $hlp 
	}
      }
    } else {

    # If tkman is already running, just update the display, else start TkMan

    set lst [winfo interps];
    if {[lsearch $lst "tkman"] == -1 } {
      set fail [catch {exec tkman &} val]

      if {$fail} {
          htmlWindow "<HTML><BODY><H1>TkMan is not available</H2></BODY></HTML>"
	  return
      }

      while {[lsearch [winfo interps] "tkman"] == -1 } {
        after 100
        }
      }
     send tkman manShowMan $manWord
  }
}



#
# And now, lets get this show on the road!
#

createDisplay

bind . <Destroy> {if {[string match "." %W]} {SaveState; exit 0}}

selectCourse [file join $Tutor(lessonHome) [format "%s.%s" $Tutor(courseName) cfg]]

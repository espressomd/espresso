###############################################################
#                                                             #
# tclline.tcl                                                 #
# ===========                                                 #
#                                                             #
# tclline: An attempt at a pure tcl readline.                 #
# the original version of this script has been taken from     #
# http://wiki.tcl.tk/16139                                    #
# and was written by Adly Abdullah (http://wiki.tcl.tk/13212) #
# who holds the Copyright of the original version             #
#                                                             #
###############################################################
# Copyright (C) 2010,2012 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#
# Use Tclx if available:
catch {
    package require Tclx

    # Prevent sigint from killing our shell:
    signal ignore SIGINT
}

# Initialise our own env variables:
foreach {var val} {
    PROMPT "> "
    HISTORY ""
    HISTORY_BUFFER 100
    COMPLETION_MATCH ""
} {
    if {![info exists env($var)]} {
        set env($var) $val
    }
}
foreach {var val} {
    CMDLINE ""
    CMDLINE_CURSOR 0
    CMDLINE_LINES 0
    HISTORY_LEVEL -1
} {
    set env($var) $val
}
unset var val

array set ALIASES {}
set forever 0

# Resource & history files:
set HISTFILE $env(HOME)/.tclline_history
set RCFILE $env(HOME)/.tcllinerc

proc ESC {} {
    return "\033"
}

proc shift {ls} {
    upvar 1 $ls LIST
    set ret [lindex $LIST 0]
    set LIST [lrange $LIST 1 end]
    return $ret
}

proc readbuf {txt} {
    upvar 1 $txt STRING

    set ret [string index $STRING 0]
    set STRING [string range $STRING 1 end]
    return $ret
}

proc goto {row {col 1}} {
    switch -- $row {
        "home" {set row 1}
    }
    print "[ESC]\[${row};${col}H" nowait
}

proc gotocol {col} {
    print "\r" nowait
    if {$col > 0} {
        print "[ESC]\[${col}C" nowait
    }
}

proc clear {} {
    print "[ESC]\[2J" nowait
    goto home
}

proc clearline {} {
    print "[ESC]\[2K\r" nowait
}

proc getColumns {} {
    set cols1 0
    set cols2 0
    if {![catch {exec stty -a} err]} {
        regexp {; *(\d+)? *columns *(\d+)?} $err -> cols1 cols2
    }
    if { [string length $cols1] > 0 } {
	return $cols1
    }
    if { [string length $cols2] > 0 } {
	return $cols2
    }
}

proc prompt {{txt ""}} {
    global env

    set prompt [subst $env(PROMPT)]
    set txt "$prompt$txt"
    foreach {end mid} $env(CMDLINE_LINES) break

    # Calculate how many extra lines we need to display.
    # Also calculate cursor position:
    set n -1
    set totalLen 0
    set cursorLen [expr {$env(CMDLINE_CURSOR)+[string length $prompt]}]
    set row 0
    set col 0
    set envcols 0

    if {$envcols <= 0} {
	set envcols 80
    }

    # Render output line-by-line to $out then copy back to $txt:
    set found 0
    set out [list]
    foreach line [split $txt "\n"] {
        set len [expr {[string length $line]+1}]
        incr totalLen $len
        if {$found == 0 && $totalLen >= $cursorLen} {
            set cursorLen [expr {$cursorLen - ($totalLen - $len)}]
            set col [expr {$cursorLen % $envcols}]
            set row [expr {$n + ($cursorLen / $envcols) + 1}]

            if {$cursorLen >= $len} {
                set col 0
                incr row
            }
            set found 1
        }
        incr n [expr {int(ceil(double($len)/$envcols))}]
        while {$len > 0} {
            lappend out [string range $line 0 [expr {$envcols-1}]]
            set line [string range $line $envcols end]
            set len [expr {$len-$envcols}]
        }
    }
    set txt [join $out "\n"]
    set row [expr {$n-$row}]

    # Reserve spaces for display:
    if {$end} {
        if {$mid} {
            print "[ESC]\[${mid}B" nowait
        }
        for {set x 0} {$x < $end} {incr x} {
            clearline
            print "[ESC]\[1A" nowait
        }
    }
    clearline
    set env(CMDLINE_LINES) $n

    # Output line(s):
    print "\r$txt"

    if {$row} {
        print "[ESC]\[${row}A" nowait
    }
    gotocol $col
    lappend env(CMDLINE_LINES) $row
}

proc print {txt {wait wait}} {
    # Sends output to stdout chunks at a time.
    # This is to prevent the terminal from
    # hanging if we output too much:
    while {[string length $txt]} {
        puts -nonewline [string range $txt 0 2047]
        set txt [string range $txt 2048 end]
        if {$wait == "wait"} {
            after 1
        }
    }
}

rename unknown _unknown
proc unknown {args} {
    global env ALIASES

    set name [lindex $args 0]
    set cmdline $env(CMDLINE)
    set cmd [string trim [regexp -inline {^\s*[^\s]+} $cmdline]]
    if {[info exists ALIASES($cmd)]} {
        set cmd [regexp -inline {^\s*[^\s]+} $ALIASES($cmd)]
    }

    set new [auto_execok $name]
    if {$new != ""} {
        set redir ""
        if {$name == $cmd && [info command $cmd] == ""} {
            set redir ">&@ stdout <@ stdin"
        }
        if {[catch {
            uplevel 1 exec $redir $new [lrange $args 1 end]} ret]
        } {
            return
        }
        return $ret
    }

    eval _unknown $args
}

proc alias {word command} {
    global ALIASES
    set ALIASES($word) $command
}

proc unalias {word} {
    global ALIASES
    array unset ALIASES $word
}

################################
# Key bindings
################################
proc handleEscapes {} {
    global env
    upvar 1 keybuffer keybuffer
    set seq ""
    set found 0
    while {[set ch [readbuf keybuffer]] != ""} {
        append seq $ch

        switch -exact -- $seq {
            "\[A" { ;# Cursor Up (cuu1,up)
                handleHistory 1
                set found 1; break
            }
            "\[B" { ;# Cursor Down
                handleHistory -1
                set found 1; break
            }
            "\[C" { ;# Cursor Right (cuf1,nd)
                if {$env(CMDLINE_CURSOR) < [string length $env(CMDLINE)]} {
                    incr env(CMDLINE_CURSOR)
                }
                set found 1; break
            }
            "\[D" { ;# Cursor Left
                if {$env(CMDLINE_CURSOR) > 0} {
                    incr env(CMDLINE_CURSOR) -1
                }
                set found 1; break
            }
            "\[H" -
            "\[7~" -
            "\[1~" { ;# home
                set env(CMDLINE_CURSOR) 0
                set found 1; break
            }
            "\[3~" { ;# delete
                if {$env(CMDLINE_CURSOR) < [string length $env(CMDLINE)]} {
                    set env(CMDLINE) [string replace $env(CMDLINE) \
                        $env(CMDLINE_CURSOR) $env(CMDLINE_CURSOR)]
                }
                set found 1; break
            }
            "\[F" -
            "\[K" -
            "\[8~" -
            "\[4~" { ;# end
                set env(CMDLINE_CURSOR) [string length $env(CMDLINE)]
                set found 1; break
            }
            "\[5~" { ;# Page Up }
            "\[6~" { ;# Page Down }
        }
    }
    return $found
}

proc handleControls {} {
  global env
  upvar 1 char char
  upvar 1 keybuffer keybuffer

  # Control chars start at a == \u0001 and count up.
  switch -exact -- $char {
      \u0001 { ;# ^a
          set env(CMDLINE_CURSOR) 0
      }
      \u0002 { ;# ^b
          if {$env(CMDLINE_CURSOR) > 0} {
      	incr env(CMDLINE_CURSOR) -1
          }
      }
      \u0004 { ;# ^d
          set env(CMDLINE) [string replace $env(CMDLINE) \
      			  $env(CMDLINE_CURSOR) $env(CMDLINE_CURSOR)]
      }
      \u0005 { ;# ^e
          set env(CMDLINE_CURSOR) [string length $env(CMDLINE)]
      }
      \u0006 { ;# ^f
          if {$env(CMDLINE_CURSOR) < [string length $env(CMDLINE)]} {
      	incr env(CMDLINE_CURSOR)
          }
      }
      \u0007 { ;# ^g
          set env(CMDLINE) ""
          set env(CMDLINE_CURSOR) 0
      }
      	\u000b { ;# ^k
          set env(YANK)  [string range $env(CMDLINE) [expr {$env(CMDLINE_CURSOR)  } ] end ]
          set env(CMDLINE) [string range $env(CMDLINE) 0 [expr {$env(CMDLINE_CURSOR) - 1 } ]]
      }
      \u0019 { ;# ^y
          if { [ info exists env(YANK) ] } {
      	set env(CMDLINE) "[string range $env(CMDLINE) 0 [expr {$env(CMDLINE_CURSOR) - 1 }]]$env(YANK)[string range $env(CMDLINE) $env(CMDLINE_CURSOR) end]"
          }
      }
      \u000e { ;# ^n
          handleHistory -1
      }
      \u0010 { ;# ^p
          handleHistory 1
      }
      \u0003 { ;# ^c
          # doExit
      }
      \u0008 -
      \u007f { ;# ^h && backspace ?
          if {$env(CMDLINE_CURSOR) > 0} {
      	incr env(CMDLINE_CURSOR) -1
      	set env(CMDLINE) [string replace $env(CMDLINE) \
      			      $env(CMDLINE_CURSOR) $env(CMDLINE_CURSOR)]
          }
      }
      \u001b { ;# ESC - handle escape sequences
          handleEscapes
      }
  }
  # Rate limiter:
  set keybuffer ""
}

proc shortMatch {maybe} {
    # Find the shortest matching substring:
    set maybe [lsort $maybe]
    set shortest [lindex $maybe 0]
    foreach x $maybe {
        while {![string match $shortest* $x]} {
            set shortest [string range $shortest 0 end-1]
        }
    }
    return $shortest
}

proc handleCompletion {} {
    global env
    set vars ""
    set cmds ""
    set execs ""
    set files ""

    # First find out what kind of word we need to complete:
    set wordstart [string last " " $env(CMDLINE) \
        [expr {$env(CMDLINE_CURSOR)-1}]]
    incr wordstart
    set wordend [string first " " $env(CMDLINE) $wordstart]
    if {$wordend == -1} {
        set wordend end
    } else {
        incr wordend -1
    }
    set word [string range $env(CMDLINE) $wordstart $wordend]

    if {[string trim $word] == ""} return

    set firstchar [string index $word 0]

    # Check if word is a variable:
    if {$firstchar == "\$"} {
        set word [string range $word 1 end]
        incr wordstart

        # Check if it is an array key:
        set x [string first "(" $word]
        if {$x != -1} {
            set v [string range $word 0 [expr {$x-1}]]
            incr x
            set word [string range $word $x end]
            incr wordstart $x
            if {[uplevel #0 "array exists $v"]} {
                set vars [uplevel #0 "array names $v $word*"]
            }
        } else {
            foreach x [uplevel #0 {info vars}] {
                if {[string match $word* $x]} {
                    lappend vars $x
                }
            }
        }
    } else {
        # Check if word is possibly a path:
        if {$firstchar == "/" || $firstchar == "." || $wordstart != 0} {
            set files [glob -nocomplain -- $word*]
        }
        if {$files == ""} {
            # Not a path then get all possibilities:
            if {$firstchar == "\[" || $wordstart == 0} {
                if {$firstchar == "\["} {
                    set word [string range $word 1 end]
                    incr wordstart
                }
                # Check executables:
                foreach dir [split $env(PATH) :] {
                    foreach f [glob -nocomplain -directory $dir -- $word*] {
                        set exe [string trimleft [string range $f \
                            [string length $dir] end] "/"]

                        if {[lsearch -exact $execs $exe] == -1} {
                            lappend execs $exe
                        }
                    }
                }
                # Check commands:
                foreach x [info commands] {
                    if {[string match $word* $x]} {
                        lappend cmds $x
                    }
                }
            } else {
                # Check commands anyway:
                foreach x [info commands] {
                    if {[string match $word* $x]} {
                        lappend cmds $x
                    }
                }
            }
        }
        if {$wordstart != 0} {
            # Check variables anyway:
            set x [string first "(" $word]
            if {$x != -1} {
                set v [string range $word 0 [expr {$x-1}]]
                incr x
                set word [string range $word $x end]
                incr wordstart $x
                if {[uplevel #0 "array exists $v"]} {
                    set vars [uplevel #0 "array names $v $word*"]
                }
            } else {
                foreach x [uplevel #0 {info vars}] {
                    if {[string match $word* $x]} {
                        lappend vars $x
                    }
                }
            }
        }
    }

    set maybe [concat $vars $cmds $execs $files]
    set shortest [shortMatch $maybe]
    if {"$word" == "$shortest"} {
        if {[llength $maybe] > 1 && $env(COMPLETION_MATCH) != $maybe} {
            set env(COMPLETION_MATCH) $maybe
            clearline
            set temp ""
            foreach {match format} {
                vars  "35"
                cmds  "1;32"
                execs "32"
                files "0"
            } {
                if {[llength [set $match]]} {
                    append temp "[ESC]\[${format}m"
                    foreach x [set $match] {
                        append temp "[file tail $x] "
                    }
                    append temp "[ESC]\[0m"
                }
            }
            print "\n$temp\n"
        }
    } else {
        if {[file isdirectory $shortest] &&
            [string index $shortest end] != "/"} {
            append shortest "/"
        }
        if {$shortest != ""} {
            set env(CMDLINE) \
                [string replace $env(CMDLINE) $wordstart $wordend $shortest]
            set env(CMDLINE_CURSOR) \
                [expr {$wordstart+[string length $shortest]}]
        } elseif {$env(COMPLETION_MATCH) != " not found "} {
            set env(COMPLETION_MATCH) " not found "
            print "\nNo match found.\n"
        }
    }
}

proc handleHistory {x} {
    global env

    set hlen [llength $env(HISTORY)]
    incr env(HISTORY_LEVEL) $x
    if {$env(HISTORY_LEVEL) > -1} {
        set env(CMDLINE) [lindex $env(HISTORY) end-$env(HISTORY_LEVEL)]
        set env(CMDLINE_CURSOR) [string length $env(CMDLINE)]
    }
    if {$env(HISTORY_LEVEL) <= -1} {
        set env(HISTORY_LEVEL) -1
        set env(CMDLINE) ""
        set env(CMDLINE_CURSOR) 0
    } elseif {$env(HISTORY_LEVEL) > $hlen} {
        set env(HISTORY_LEVEL) $hlen
    }
}

################################
# History handling functions
################################

proc getHistory {} {
    global env
    return $env(HISTORY)
}

proc setHistory {hlist} {
    global env
    set env(HISTORY) $hlist
}

proc appendHistory {cmdline} {
    global env
    set old [lsearch -exact $env(HISTORY) $cmdline]
    if {$old != -1} {
        set env(HISTORY) [lreplace $env(HISTORY) $old $old]
    }
    lappend env(HISTORY) $cmdline
    set env(HISTORY) \
        [lrange $env(HISTORY) end-$env(HISTORY_BUFFER) end]
}

################################
# main()
################################

proc rawInput {} {
    fconfigure stdin -buffering none -blocking 0
    fconfigure stdout -buffering none -translation crlf
    exec stty raw -echo
}

proc lineInput {} {
    fconfigure stdin -buffering line -blocking 1
    fconfigure stdout -buffering line
    exec stty -raw echo
}

proc doExit {{code 0}} {
    global env HISTFILE

    # Reset terminal:
    print "[ESC]c[ESC]\[2J" nowait
    lineInput

    set hlist [getHistory]
    if {[llength $hlist] > 0} {
        set f [open $HISTFILE w]
        foreach x $hlist {
            # Escape newlines:
            puts $f [string map {
                \n "\\n"
                "\\" "\\b"
            } $x]
        }
        close $f
    }

    exit $code
}

if {[file exists $RCFILE]} {
    source $RCFILE
}

# Load history if available:
if {[llength $env(HISTORY)] == 0} {
    if {[file exists $HISTFILE]} {
        set f [open $HISTFILE r]
        set hlist [list]
        foreach x [split [read $f] "\n"] {
            if {$x != ""} {
                # Undo newline escapes:
                lappend hlist [string map {
                    "\\n" \n
                    "\\\\" "\\"
                    "\\b" "\\"
                } $x]
            }
        }
        setHistory $hlist
        unset hlist
        close $f
    }
}

rawInput

# This is to restore the environment on exit:
# Do not unalias this!
alias exit doExit

proc tclline {} {
    global env
    set char ""
    set keybuffer [read stdin]
    set env(COLUMNS) [getColumns]

    while {$keybuffer != ""} {
        if {[eof stdin]} return
        set char [readbuf keybuffer]
        if {$char == ""} {
            # Sleep for a bit to reduce CPU time:
            after 40
            continue
        }

        if {[string is print $char]} {
            set x $env(CMDLINE_CURSOR)

            if {$x < 1 && [string trim $char] == ""} continue

            set trailing [string range $env(CMDLINE) $x end]
            set env(CMDLINE) [string replace $env(CMDLINE) $x end]
            append env(CMDLINE) $char
            append env(CMDLINE) $trailing
            incr env(CMDLINE_CURSOR)
        } elseif {$char == "\t"} {
            handleCompletion
        } elseif {$char == "\n" || $char == "\r"} {
            if {[info complete $env(CMDLINE)] &&
                [string index $env(CMDLINE) end] != "\\"} {
                lineInput
                print "\n" nowait
                uplevel #0 {
                    global env ALIASES

                    # Handle aliases:
                    set cmdline $env(CMDLINE)
                    set cmd [string trim [regexp -inline {^\s*[^\s]+} $cmdline]]
                    if {[info exists ALIASES($cmd)]} {
                        regsub -- "(?q)$cmd" $cmdline $ALIASES($cmd) cmdline
                    }

                    # Perform glob substitutions:
                    set cmdline [string map {
                        "\\*" \0
                        "\\~" \1
                    } $cmdline]
                    while {[regexp -indices \
                        {([\w/\.]*(?:~|\*)[\w/\.]*)+} $cmdline x]
                    } {
                        foreach {i n} $x break
                        set s [string range $cmdline $i $n]
                        set x [glob -nocomplain -- $s]

                        # If glob can't find anything then don't do
                        # glob substitution, pass * or ~ as literals:
                        if {$x == ""} {
                            set x [string map {
                                "*" \0
                                "~" \1
                            } $s]
                        }
                        set cmdline [string replace $cmdline $i $n $x]
                    }
                    set cmdline [string map {
                        \0 "*"
                        \1 "~"
                    } $cmdline]

                    # Run the command:
                    catch $cmdline res
                    if {$res != ""} {
                        print "$res\n"
                    }

                    # Append HISTORY:
                    set env(HISTORY_LEVEL) -1
                    appendHistory $env(CMDLINE)

                    set env(CMDLINE) ""
                    set env(CMDLINE_CURSOR) 0
                    set env(CMDLINE_LINES) {0 0}
                }
                rawInput
            } else {
                set x $env(CMDLINE_CURSOR)

                if {$x < 1 && [string trim $char] == ""} continue

                set trailing [string range $env(CMDLINE) $x end]
                set env(CMDLINE) [string replace $env(CMDLINE) $x end]
                append env(CMDLINE) $char
                append env(CMDLINE) $trailing
                incr env(CMDLINE_CURSOR)
            }
        } else {
            handleControls
        }
    }
    prompt $env(CMDLINE)
}
tclline

fileevent stdin readable tclline
vwait forever
doExit

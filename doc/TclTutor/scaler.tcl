##############################################################
# proc Scaler {parent row}--
#    Create a toolbutton icon that will evaluate a command when 
#    clicked and display an information label if the cursor 
#    rests on it
#
# Arguments
#   parent	The parent window for this scaler
#   row		The row to scale.
#
# Results
#   Creates a new Scaler and returns the name to the 
#   calling script
#   Throws an error if unrecognized arguments or required args
#   are missing.
#

proc makeScaler {parent row background {text {}}} {
    global Scalers
    
    set info "WindowSizer - Hold down button one and move mouse up to grow, down to shrink"
    
    if {![info exists Scalers(unique)]} {set Scalers(unique) 0}
    set name $parent-s$Scalers(unique)
    incr Scalers(unique)
   
    canvas $name -height 12 -width 400 -background $background \
	    -highlightbackground $background -relief flat
    $name create polygon 0 5 5 10 10 5 5 0 -fill blue -tag mover
    $name create text 40 6 -anchor w -text $text
#    focus $name

    $name bind mover <B1-Motion> "changeScale $name $row %Y"
    

    # Bind the background to change color when the 
    #    cursor enters

    $name bind mover <Enter> \
        "$name itemconfigure mover -fill green"

    # Bind the background to change color when the cursor 
    #  leaves, and cancel any pending display of an info label.

    $name bind mover <Leave> \
        "$name itemconfigure mover -fill blue; \
       [namespace current]::toolInfo::cancelInfo \
       $name;"

    # Whenever the button moves, reset the Motion task.

    $name bind  mover <Motion> \
        "[list [namespace current]::toolInfo::scheduleInfo \
        $name $info ]"

    # And return the name

    return $name
}

proc changeScale {id row y } {
    global Scalers

    if {[info exists Scalers($id.y)]} {
        set v [expr $y -$Scalers($id.y)]
	if {$v > 0} {set v -1} else {set v 1}
	incr Scalers($row) $v
	if {$Scalers($row) < 0} {set Scalers($row) 0}

	grid rowconfigure . $row -weight $Scalers($row)
    }
    set Scalers($id.y) $y
}

#
# The infoIdentifier, scheduleInfo, cancelInfo, createInfo 
#  and deleteInfo procedures are created in the toolInfo
#  namespace to
#  1) avoid polluting the global space.
#  2) make the  "afters" array of identifiers returned by the
#     after command private to these procedures.
#  3) make the "afters" array persistent
#  4) make the  "binding" array of identifiers returned by 
#     the after command private to these procedures.
#  5) make the "binding" array persistent
#

namespace eval toolInfo {
    
##############################################################
    # proc infoIdentifier {toolButtonName}--
    #    Convert a button path name into a string that can 
    #    be used to identify the label.
    #    The label will be placed in the root window,
    #    independent of frames that may hold the toolButtons.
    #    This allows the labels to overlap other widgets 
    #    without causing the frame to resize itself.
    #    Because of this, the label can only have one "."
    #    in its name, and other "." are converted to "_"
    #
    # Arguments
    #    toolButtonNameThe full window name of a toolButton
    #
    # Results
    #   Returns the modified string

    proc infoIdentifier {toolButtonName} {
        regsub -all {\.} $toolButtonName "_" infoID
        return $infoID
    }

    
##############################################################
    # proc cancelInfo {toolButtonName}--
    #    Cancels the scheduled display of an Info label
    # Arguments
    #   toolButtonNameThe name of the toolbutton
    #
    # Results
    #   If this label was scheduled to be displayed, 

    #    it is disabled.

    proc cancelInfo {toolButtonName} {
        variable afters
        if {[info exists afters($toolButtonName)]} {
            after cancel $afters($toolButtonName);
        }
    }

    
##############################################################
    #     proc scheduleInfo {name info }--
    #    Cancel any existing "after" scripts for this button.
    #
    # Arguments
    #   nameThe name of this toolbutton
    #   infoThe text to display in the  info label
    #
    # Results
    #   Places an item in the "after" queue to be evaluated
    #      2 seconds from now

    proc scheduleInfo {toolButtonName info } {
        variable afters

        cancelInfo $toolButtonName

        set afters($toolButtonName) \
            [after 1500 [list [namespace current]::createInfo \
            $toolButtonName $info ]]
    }

    
##############################################################
    #     proc createInfo { toolButtonName info}--
    #    Creates an info label, and sets up bindings 
    #    to delete it.
    #
    # Arguments
    #   toolButtonName  The name of the parent button
    #   info            The information to display.
    #
    # Results
    #   Creates and places a label.
    #   Adds a script to the motion binding for the parent 
    #   button,
    #   and clears it when the motion occurs.

    proc createInfo { toolButtonName info} {

        variable binding

        # Create a name for this label, and create the label

        set labelName ".toolbutton_[infoIdentifier \
            $toolButtonName]"
        label $labelName -text $info

        # Determine the location for this label.
        #
        # winfo returns the geometry of the window within 
        #   its parent as WIDTHxHEIGHT+XOFFSET+YOFFSET
        #   split at the + to separate out the X and Y 
        #   offsets.

        set coordLst [split [winfo geometry \
            $toolButtonName] "+"]
        set x [lindex $coordLst 1]
        set y [lindex $coordLst 2]

        # Place the label over the appropriate toolButton

        place $labelName -x [expr $x+6] -y [expr $y+10]

        # If this toolButton is close to the right margin 
        #   of the root window, the label won't be
        #   completely displayed.
        # In this case, the label must be moved to the 
        #   left to fit on the window.

        # Update to make the label appear.  The $labelName 
        #  geometry is not valid until this update occurs.

        update idle

        # Get the width of the root window, and the width 
        #   of the $labelName label.
        # If the X location + width goes over the edge,
        # move the label to the left until it fits.

        set rootWidth [lindex [split [winfo geometry .] "x"] 0]
        set labelWidth [lindex [split \
            [winfo geometry $labelName] "x"] 0]

        if {[expr $x + $labelWidth] > $rootWidth} {
            set x [expr $rootWidth - $labelWidth]
            place $labelName -x $x -y [expr $y+10]
        }

        # Save the old event binding, and set a new one,
        # that binds a Motion event to destroying the label

        set binding($toolButtonName) \
            [bind $toolButtonName <Motion>]

        bind $toolButtonName <Motion> \
            "[namespace current]::deleteInfo $toolButtonName";

        bind $labelName <Motion> \
            "[namespace current]::deleteInfo $toolButtonName";
    }

    
##############################################################
    #     proc deleteInfo {toolButtonName}--
    #    Deletes the label associated with this button, and
    #    resets the binding to the original actions.
    # Arguments
    #   toolButtonName    The name of the parent button
    #
    # Results
    #   No more label, and bindings as they were before the
    #   mouse paused
    #   on this button

    proc deleteInfo {toolButtonName} {
        variable binding

        set labelName \
            ".toolbutton_[infoIdentifier $toolButtonName]"
        destroy $labelName;

        bind $toolButtonName <Motion> \
            $binding($toolButtonName);
    }
}

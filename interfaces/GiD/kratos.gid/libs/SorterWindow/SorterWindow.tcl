namespace eval SorterWindow {
    variable winpath
    variable plot
    variable PosX
    variable PosY
    variable Delta
    variable Index
    variable Listado
    variable Scrollposition
    variable Scrollarea
    variable imgs
    variable window_state
    variable data_source
    variable data_source_list
}

proc SorterWindow::Init {} {
    variable imgs
    variable window_state
    variable winpath
    set winpath ".gid.sortwindow"
    #variable dir
    set dir [file dirname [info script]]
    set imgs(lock) [image create photo -file [file join $dir images lock.png]]
    set imgs(unlock) [image create photo -file [file join $dir images unlock.png]]
    set window_state "unlocked"
    
    variable data_source
    set data_source "Origin"
    variable data_source_list
    set data_source_list [list "Origin" "Source"]
}

#
# Open a widow for the test dialog. Left half = Scrollbox.
# Right: Some entries and a message for testing the code.
#
proc SorterWindow::Window { } {
    variable PosX
    variable PosY
    variable Delta
    variable winpath
    set w $winpath
    
    catch {destroy $w}
    
    ###################
    # CREATING WIDGETS
    ###################
    toplevel $w -class Toplevel -relief groove 
    wm maxsize $w 500 300
    wm minsize $w 500 300
    wm overrideredirect $w 0
    wm resizable $w 1 1
    wm deiconify $w
    wm title $w [= "Sort conditions window"]
    set Delta 30
    
    SorterWindow::RefreshWindow
    
}

proc SorterWindow::RefreshWindow { } {
    variable winpath
    set w $winpath
    if {[winfo exists $w.fr1]} {destroy $w.fr1}
    if {[winfo exists $w.fr2]} {destroy $w.fr2}
    if {[winfo exists $w.buts]} {destroy $w.buts}
    set fr1 [ttk::frame $w.fr1]
    set fr2 [ttk::labelframe $w.fr2 -text [= "Sort the conditions:"]]
    set buts [ttk::frame $w.buts -style BottomFrame.TFrame]
    
    # One call creates the listbox
    #
    set Canv [SortByDragListbox $fr1 20 20 150 240]
    #
    # Widgets on the right side
    #
    ttk::button $buts.q -text Cancel -command [list destroy $w] -style BottomFrame.TButton
    ttk::button $buts.ok -text Ok -command [list destroy $w] -style BottomFrame.TButton
    SorterWindow::ConfigurationFrame $fr2
    
    
    grid $buts.ok $buts.q -sticky sew
    
    grid $fr1 -sticky nsew -row 0 -column 0
    grid $fr2 -sticky nw -row 0 -column 1 -padx 20
    grid $buts -sticky sew -columnspan 2
    grid columnconfigure $w 0 -weight 1 
    grid rowconfigure $w 0 -weight 1
    if { $::tcl_version >= 8.5 } { grid anchor $buts center }
}

proc SorterWindow::ConfigurationFrame {w} {
    variable imgs
    variable window_state
    variable data_source_list
    
    set locktext [= "Drag&Drop unlocked"]:
    set unlocktext [= "Drag&Drop locked"]:
    set lt $locktext
    if {$window_state eq "locked"} {set lt $unlocktext}
    set lab1 [ttk::label $w.l1 -text $lt]
    
    set li $imgs(unlock)
    if {$window_state eq "locked"} {set li $imgs(lock)}
    set b1 [ttk::button $w.b1 -image $li -command [list SorterWindow::LockButtonClicked]]
    grid $lab1 $b1 -sticky w
    
    set lab2 [ttk::label $w.l2 -text [= "Data source"]:]
    set cb2 [ttk::combobox $w.cb2 -textvariable SorterWindow::data_source -values $data_source_list -width 10 -state readonly]
    bind $cb2 <<ComboboxSelected>> SorterWindow::DataSourceChanged
    grid $lab2 $cb2 -sticky w
}

proc SorterWindow::LockButtonClicked { } {
    variable window_state
    if {$window_state eq "locked"} {set window_state "unlocked"} {set window_state "locked"}
    SorterWindow::RefreshWindow
}
proc SorterWindow::DataSourceChanged { } {
    SorterWindow::RefreshWindow
}

#
# Create a pseudo-listbox with canvas elements. Looks like a listbox,
# but is really a canvas, and all widgets only pretend to be what they seem.
#
proc SorterWindow::SortByDragListbox { w XNull YNull width height } {
    variable Index
    variable Listado
    variable Scrollposition
    variable Scrollarea
    variable window_state
    catch {destroy $w.cv}
    set Canv [canvas $w.cv -borderwidth 0 -highlightthickness 0  -height [expr $height + 2*$YNull] -width [expr $width + 2*$XNull] ]

    $Canv create rectangle $XNull $YNull [expr $XNull + $width] [expr $YNull + $height] -outline black -width 1 -fill white -tags Box
    
    if {$window_state eq "locked"} {
        $Canv create rectangle [expr $XNull - 1] [expr $YNull - 1] [expr $XNull + $width + 1] [expr $YNull + $height + 1] -outline red -width 1 -tags Box -fill seashell
        $Canv configure -state disabled
    } else {
        $Canv create rectangle [expr $XNull - 1] [expr $YNull - 1] [expr $XNull + $width + 1] [expr $YNull + $height + 1] -outline grey50 -width 1 -tags Box
    }
    
    catch {destroy $w.lbscroll}
    ttk::scrollbar $w.lbscroll -command "SorterWindow::SortByDragListboxScroll $w" -orient vert
    grid $Canv -row 0 -column 0 -sticky nw
    grid $w.lbscroll -row 0 -column 1 -sticky wns
    
    SorterWindow::FillListado
    set Scrollposition 0
    set Schrifthoehe 16
    # Cambiar el 20 por el numero de tal
    set Scrollarea [expr $Schrifthoehe * 20]
    SorterWindow::SortByDragListboxScroll $w scroll 0.0 units
    
    return $Canv
}

#
# Scrollbar code.
#
proc SorterWindow::SortByDragListboxScroll { w {was moveto} {Zahl 0.0} {Einheit units} } {
    variable Index
    variable Listado
    variable Scrollposition
    variable Scrollarea
    
    set Canv  $w.cv
    set height [lindex [$Canv configure -height] 4]
    set Schrifthoehe 16
    set Scrollposition 0
    
    if {$was == "scroll"} {
        if {$Einheit == "pages"} {
            incr Scrollposition [expr int($Zahl * $height - 20)]
        } else {
            incr Scrollposition [expr 20 * int($Zahl)]
        }
    } else {
        set Scrollposition [expr int($Zahl * $Scrollarea)]
    }
    
    # Limit the scrollposition to sensible values
    #
    if {$Scrollposition > [expr $Scrollarea - $height]} {
        set Scrollposition [expr $Scrollarea - $height]
    }
    if {$Scrollposition < 0} {set Scrollposition 0}
    #
    # Delete Index and built anew from scratch. In priciple all entries could
    # be moved, but this is messy at the edges.
    #
    set yPos [expr 32 - $Scrollposition]
    for {set i 0} {$i < [array size Listado]} {incr i} {
        $Canv delete ent$i
        if {$yPos < 20} {
            incr yPos $Schrifthoehe
            continue
        }
        #
        if {$yPos < [expr $height - 20]} {
            $Canv create text 24 $yPos -text $Listado($Index($i)) -anchor w -fill black -tags ent$i
            incr yPos $Schrifthoehe
            
            $Canv bind ent$i <1>     "SorterWindow::plotDown $Canv %x %y"
            $Canv bind ent$i <B1-Motion>    "SorterWindow::plotMove $Canv %x %y"
            $Canv bind ent$i <ButtonRelease-1> "SorterWindow::plotCopy $w $Canv %x %y $i"
        }
    }
    #
    $w.lbscroll set [expr double($Scrollposition) / $Scrollarea] [expr double($height + $Scrollposition) / $Scrollarea]
}

#
# plotDown --
# This procedure is invoked when the mouse is pressed over one of the
# data points. It sets up state to allow the point to be dragged.
#
# Arguments:
# w -       The canvas window.
# x, y -    The coordinates of the mouse press.
#
proc SorterWindow::plotDown {w x y} {
    variable plot
    #
    $w dtag selected
    $w addtag selected withtag current
    $w raise current
    set plot(lastX) $x
    set plot(lastY) $y
}

# plotMove --
# This procedure is invoked during mouse motion events. It drags the
# current item.
#
# Arguments:
# w -       The canvas window.
# x, y -    The coordinates of the mouse.
#
proc SorterWindow::plotMove { w x y } {
    variable plot
    variable PosX
    variable PosY
    
    $w move selected [expr $x-$plot(lastX)] [expr $y-$plot(lastY)]
    set plot(lastX) $x
    set plot(lastY) $y
    set PosX        $x
    set PosY        $y
}

#
# When the mouse button is released, this routine determines the new
# position and re-orders the list.
#
proc SorterWindow::plotCopy { w Cv x y i } {
    variable Delta
    variable Index
    variable Scrollposition
    variable Scrollarea
    
    set Schrifthoehe 16
    #
    # Get the new position. Delta is a fudge factor which is different
    # between different operating systems.
    #
    set Rang [expr int(($y - $Delta + $Scrollposition) / $Schrifthoehe)]
    set Temp $Index($i)
    if {$Rang > $i} {
        for {set j $i} {$j < $Rang} {incr j} {
            set Index($j) $Index([expr $j + 1])
        }
    } elseif {$Rang == $i} {
        set number [expr double($Scrollposition) / $Scrollarea]
        SorterWindow::SortByDragListboxScroll $w scroll $number units
        return
    } else {
        set Rang [expr $Rang + 1]
        for {set j $i} {$j > $Rang} {incr j -1} {
            set Index($j) $Index([expr $j - 1])
        }
    }
    set Index($Rang) $Temp
    #
    # Now scroll the list to the right position.
    #
    set number [expr double($Scrollposition) / $Scrollarea]
    SorterWindow::SortByDragListboxScroll $w scroll $number units
}

proc SorterWindow::FillListado { } {
    variable Index
    variable Listado
    variable data_source
    catch {unset Index }
    set Items [list ]
    
    set word "Item"
    if {$data_source eq "Origin"} {set word "Ejemplar"}
    for {set i 1} {$i < 21} {incr i} {
        lappend Items "$word $i"
    }
    for {set i 0} {$i < [llength $Items]} {incr i} {
        set Listado($i) "[lindex $Items $i]"        
        lappend Index($i) $i        
    }
    
}


proc SorterWindow::SorterWindow {{init default}} {
    
    if {[GiD_Info problemtypepath] ne "E:/PROYECTOS/Kratos/interfaces/GiD/kratos.gid"} {
        catch {SorterWindow::TMP_Cards $w}
    } {
        #SorterWindow::Window show $w
        SorterWindow::Window 
    }
}

proc SorterWindow::TMP_Cards { } {
    catch {destroy .gid.w1.c}
    catch {destroy .gid.w1}
    toplevel .gid.w1
    pack [canvas .gid.w1.c -bg darkgreen -borderwidth 0 -highlightthickness 0] -expand 1 -fill both
    
    image create photo ::SorterWindow::ah -format gif -data {
        R0lGODlhRwBgAKEAAH//1AAAAP////8AACH5BAEAAAAALAAAAABHAGAAAAL/BIKpy+0PYzBH2I
        uz3rz7L0wBSJamiZzquqbawMZyOL7zfbrYwOP+p7vwYL+iJijo9YxMWkZJbBaDy2RUiqMOh9if
        bgvuZmvW51XMcm2FXHRM3bZW3SokfUq+G+36MRsW15dTs7YmOGgBZnhYAqd4xujhqBjZSPZYab
        mzmAmUJ9ep+RQqStryaVqaioK66umKCKsqK9lKe2R7i8Gnu5vb6wTMwStMLLPI6WdESen1u4KJ
        6WMM/Sit/GN9fUOtot2M7fMdNv1cPY7HND7HbX5uvef+Tp4ute3cRR8vFrjP39XtVkBaA2UVhH
        XQVcJVC1M1NPWQVMRQEztVzHTxDqSMYm76ceTH6SOWayKbaLNQUh28YLROspRVqE1KlYCqzLxz
        k05ONztnIJMp79DPJT19nrEZVOjKl7C2FTXKxlcvKBmeQmVn1ejGpJH6MW25IavPsFwZlnVYQV
        hVChLaun3b1kABADs=}
    image create photo ::SorterWindow::as -format gif -data {
        R0lGODlhRwBgAKEAAH//1AAAAP///////yH5BAEAAAAALAAAAABHAGAAAAL/BIKpy+0PYzBH2I
        uz3rz7L0wBSJamiZzquqbZyMZyCGP1jJfuleQ+uLP0fsRNMBUsFo+jpPK3i96ePumCuoQ9sFDt
        1MllIYc0cPg0tk7PKjO7un535VRnnI7+uvHA25XfUtMA2OY1SKhjyICYKOSQcwfH00QDuTeTRI
        apUKcXibKodBnDyZn1VFr2GSja46qJM5qXRjvXRfsq23fblNur6wEc7BuiarxpWUtsSrrap9xr
        zMwqU/paDC1s49xhPVZW6c2toT2ZDb4MS1Iu7Ut5JV4YiPuLq1oLOn+fi02fllfoHzx099Cp6z
        bOXMFvDKPBA6bNnTSC4n7lA1XPkUCD3y68NcLI0d5DfgeNJOTYLkJKcOtOkhQJ02JJci7xLcRW
        UKPFlrP2jQw5aOaLmvUWbYHAchfGou6MBhW6LaBORSGnary41FPHrbD+8AyI9OrUR1hnkc3pc9
        pXMaHMZWoLsJNXIuxazrXV6hDdk4m27GVkCXAsvoKHFrZr+DAQmooR9pvU2OShunL8Un5jmTAf
        vZcrT+vsWU/kYIxHlzX9ATQj1XIFDWGNiowp2LETV0Kd1jXusX4005ESdXc40cK/BS9uxzcedb
        S5xGmO5bnywtB/VxDOYYIBCdy7e+duoAAAOw==}
    
    bind .gid.w1.c <ButtonPress-1> { SorterWindow::bPress1 %W %x %y}
    bind .gid.w1 <Configure> {wm title .gid.w1 "Merry Christmas"}
}
proc  SorterWindow::bPress1 {w x y} {
    
    set i [lindex [$w find overlapping $x $y $x $y] end]
    if { $i == "" } {
        set i [$w create image $x $y -anchor nw -image ::SorterWindow::ah -tags card]
        SorterWindow::setFocus $w $i
    } else {
        SorterWindow::setFocus $w $i
    }
    
}


proc SorterWindow::setFocus {w i} {
    
    if { [.gid.w1.c find withtag hasfocus] != "" } {
        $w itemconfigure [.gid.w1.c find withtag hasfocus] -image ::SorterWindow::ah
        $w dtag [.gid.w1.c find withtag hasfocus] hasfocus
    }
    $w focus $i
    $w itemconfigure $i -image ::SorterWindow::as
    $w addtag hasfocus withtag $i
    
}

SorterWindow::Init
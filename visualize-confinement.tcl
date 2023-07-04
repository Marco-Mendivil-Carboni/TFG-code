set r_LJ 5.0
set res 32

if {$argc==3} {
    package require topotools
    set sim_dir [lindex $argv 0]
    set sim_idx [lindex $argv 1]
    set file1 [format "%s/initial-configuration-%03d.gro" $sim_dir $sim_idx]
    set file2 [format "%s/trajectory-positions-%03d.trr" $sim_dir $sim_idx]
    color Display Background white
    display cuedensity 0.2
    mol new $file1 autobonds off
    set N [molinfo top get numatoms]
    for {set i 0} {$i<[expr $N-1]} {incr i} {
        topo addbond $i [expr $i+1]
    }
    set sel [atomselect top "type X"]
    $sel set radius $r_LJ
    mol modstyle 0 top VDW 1.0 $res
    # mol modstyle 0 top licorice $r_LJ $res $res
    # mol modstyle 0 top CPK 4.0 [expr 4.0*$r_LJ/2.0] $res $res
    mol addfile $file2
} else {
    puts "You forgot the input."
    exit
}

draw material Transparent
set R [lindex $argv 2]
draw sphere {0 0 0} radius $R resolution $res

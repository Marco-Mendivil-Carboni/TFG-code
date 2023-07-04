set ::pi 3.1415926535897931
set r_LJ 5.0
set r_p 11.5
set l_p 11.9
set res 32

if {$argc==2} {
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

set ::ex {1.0 0.0 0.0}
set ::ey {0.0 1.0 0.0}
set ::ez {0.0 0.0 1.0}
proc draw_ring { res x_ctr_a x_ctr_b r_ctr_a r_ctr_b theta_a theta_b side} {
    for {set j 0} {$j < $res} {incr j} {
        set theta_c [expr 2.0*$::pi*($j+0.0)/$res]
        set theta_d [expr 2.0*$::pi*($j+1.0)/$res]
        set x_v_0 [vecscale $::ex [expr $side*$x_ctr_a]]
        set x_v_1 [vecscale $::ex [expr $side*$x_ctr_a]]
        set x_v_2 [vecscale $::ex [expr $side*$x_ctr_b]]
        set x_v_3 [vecscale $::ex [expr $side*$x_ctr_b]]
        set y_v_0 [vecscale $::ey [expr $r_ctr_a*cos($theta_c)]]
        set y_v_1 [vecscale $::ey [expr $r_ctr_a*cos($theta_d)]]
        set y_v_2 [vecscale $::ey [expr $r_ctr_b*cos($theta_c)]]
        set y_v_3 [vecscale $::ey [expr $r_ctr_b*cos($theta_d)]]
        set z_v_0 [vecscale $::ez [expr $r_ctr_a*sin($theta_c)]]
        set z_v_1 [vecscale $::ez [expr $r_ctr_a*sin($theta_d)]]
        set z_v_2 [vecscale $::ez [expr $r_ctr_b*sin($theta_c)]]
        set z_v_3 [vecscale $::ez [expr $r_ctr_b*sin($theta_d)]]
        set v_0 [vecadd $x_v_0 $y_v_0 $z_v_0]
        set v_1 [vecadd $x_v_1 $y_v_1 $z_v_1]
        set v_2 [vecadd $x_v_2 $y_v_2 $z_v_2]
        set v_3 [vecadd $x_v_3 $y_v_3 $z_v_3]
        set x_n_0 [vecscale $::ex [expr -$side*sin($theta_a)]]
        set x_n_1 [vecscale $::ex [expr -$side*sin($theta_a)]]
        set x_n_2 [vecscale $::ex [expr -$side*sin($theta_b)]]
        set x_n_3 [vecscale $::ex [expr -$side*sin($theta_b)]]
        set y_n_0 [vecscale $::ey [expr cos($theta_c)*cos($theta_a)]]
        set y_n_1 [vecscale $::ey [expr cos($theta_d)*cos($theta_a)]]
        set y_n_2 [vecscale $::ey [expr cos($theta_c)*cos($theta_b)]]
        set y_n_3 [vecscale $::ey [expr cos($theta_d)*cos($theta_b)]]
        set z_n_0 [vecscale $::ez [expr sin($theta_c)*cos($theta_a)]]
        set z_n_1 [vecscale $::ez [expr sin($theta_d)*cos($theta_a)]]
        set z_n_2 [vecscale $::ez [expr sin($theta_c)*cos($theta_b)]]
        set z_n_3 [vecscale $::ez [expr sin($theta_d)*cos($theta_b)]]
        set n_0 [vecadd $x_n_0 $y_n_0 $z_n_0]
        set n_1 [vecadd $x_n_1 $y_n_1 $z_n_1]
        set n_2 [vecadd $x_n_2 $y_n_2 $z_n_2]
        set n_3 [vecadd $x_n_3 $y_n_3 $z_n_3]
        if { $side > 0 } {
            draw trinorm $v_0 $v_1 $v_2 $n_0 $n_1 $n_2
            draw trinorm $v_3 $v_2 $v_1 $n_3 $n_2 $n_1
        } else {
            draw trinorm $v_2 $v_1 $v_0 $n_2 $n_1 $n_0
            draw trinorm $v_1 $v_2 $v_3 $n_1 $n_2 $n_3
        }
    }
}

set r_cyl [expr $r_p-$r_LJ]
draw_ring $res 0.0 [expr $l_p/2.0] $r_cyl $r_cyl 0.0 0.0 +1
draw_ring $res 0.0 [expr $l_p/2.0] $r_cyl $r_cyl 0.0 0.0 -1
for {set i 0} {$i < [expr $res/4.0]} {incr i} {
    set theta_a [expr 2.0*$::pi*($i+0.0)/$res]
    set theta_b [expr 2.0*$::pi*($i+1.0)/$res]
    set x_ctr_a [expr $l_p/2.0+$r_LJ*sin($theta_a)]
    set x_ctr_b [expr $l_p/2.0+$r_LJ*sin($theta_b)]
    set r_ctr_a [expr $r_p-$r_LJ*cos($theta_a)]
    set r_ctr_b [expr $r_p-$r_LJ*cos($theta_b)]
    draw_ring $res $x_ctr_a $x_ctr_b $r_ctr_a $r_ctr_b $theta_a $theta_b +1
    draw_ring $res $x_ctr_a $x_ctr_b $r_ctr_a $r_ctr_b $theta_a $theta_b -1
}
set theta_a $theta_b
set x_ctr_a $x_ctr_b
set r_ctr_a $r_ctr_b
set r_ctr_b [expr 2.0*$r_p]
draw_ring $res $x_ctr_a $x_ctr_b $r_ctr_a $r_ctr_b $theta_a $theta_b +1
draw_ring $res $x_ctr_a $x_ctr_b $r_ctr_a $r_ctr_b $theta_a $theta_b -1

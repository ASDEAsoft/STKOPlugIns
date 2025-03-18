
# a list of all monitor and custom function actors to be called by the MonitorFunction
set all_custom_functions {}
set all_monitor_actors {}

# the main custom function caller that will call all actors in $all_monitor_actors and in $all_custom_functions list
proc CustomFunctionCaller {step_id dt T n_iter norm perc process_id is_parallel} {
	global all_monitor_actors
	global all_custom_functions
	# Call monitors: we pass the parameters needed
	foreach p $all_monitor_actors {
		$p $step_id $dt $T $n_iter $norm $perc $process_id $is_parallel
	}
	# Call all other custom functions
	foreach p $all_custom_functions {
		$p
	}
}

# Misc_commands region

region 1 \
-nodeRange 1 8

region 2 \
-node 5

region 3 \
-nodeRange 1 8

region 4 \
-ele 4

region 5 \
-ele 7

region 6 \
-ele 5

region 7 \
-ele 6

#TCL script: Recorder (7)
recorder Node -file ./out/1000C_2497kGravity_Cyclic/node_disp.out -region 1 -dof 1 disp

#TCL script: Recorder (7)
recorder Node -file ./out/1000C_2497kGravity_Cyclic/col_force.out -region 2 -dof 1 reaction

#TCL script: Recorder (7)
recorder Node -file ./out/1000C_2497kGravity_Cyclic/rotation_top.out -region 3 -dof 5 disp

#TCL script: Recorder (7)
recorder Element -file ./out/1000C_2497kGravity_Cyclic/fail.out -time -region 4 material 6 failFlag

#TCL script: Recorder (7)
recorder Element -file ./out/1000C_2497kGravity_Cyclic/force_fiber.out -region 5 force

#TCL script: Recorder (7)
recorder Element -file ./out/1000C_2497kGravity_Cyclic/deformation_fiber.out -region 6 deformation

#TCL script: Recorder (7)
recorder Element -file ./out/1000C_2497kGravity_Cyclic/rotation_j.out -region 7 material 6 rotationAndLimits

# Constraints.sp fix
	fix 5 1 1 1 1 1 1
	fix 6 0 1 0 0 0 1

# Patterns.addPattern loadPattern
pattern Plain 16 1 {

# Loads.Force NodeForce
	load 6 0.0 0.0 -2497.0 0.0 0.0 0.0
}

# Patterns.addPattern loadPattern
pattern Plain 17 1 {

# Loads.Force NodeForce
	load 6 1.0 0.0 0.0 0.0 0.0 0.0
}

# Done!
puts "ANALYSIS SUCCESSFULLY FINISHED"

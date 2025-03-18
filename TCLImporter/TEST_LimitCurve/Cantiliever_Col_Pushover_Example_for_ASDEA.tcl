### Title:  1000C_2497kGravity_Cyclic
### Author: A. D. Sen, P. Ranchal (Modified 10/7/2021)
wipe;

puts "Starting Analysis of Model: 1000C_2497kGravity_Cyclic"

####################################################################################
## Model
####################################################################################
model BasicBuilder -ndm 3 -ndf 6;

############# Nodes #############
node 1 0.00000 0.00000 0.00000;
node 2 0.00000 0.00000 75.00000;
node 3 0.00000 0.00000 0.00000;
node 4 0.00000 0.00000 0.00000;
node 5 0.00000 0.00000 8.0;
node 6 0.00000 0.00000 75.00000;
node 7 0.00000 0.00000 75.00000;
node 8 0.00000 0.00000 67;
fix 1 1 1 1 1 1 1;
fix 2 0 1 0 0 0 1;

############# Materials #############
# uniaxialMaterial 	Concrete01 $matTag 	$fpc 	$epsc0 		$fpcu 		$epsU
uniaxialMaterial 	Concrete01 1 		-7.50 	-0.02511 	-7.4925 	-0.14630;
uniaxialMaterial 	Concrete01 4 		-7.50 	-0.02511 	-7.4925 	-0.14630;

# uniaxialMaterial 	Concrete02 $matTag 	$fpc 	$epsc0 		$fpcu 		$epsU 		$lambda 	$ft 		$Ets
uniaxialMaterial 	Concrete02 2 		-7.50 	-0.02511 	-7.4925 	-0.14630 	0.10000 	0.27698 	38.18710;
uniaxialMaterial 	Concrete02 5 		-7.50 	-0.02511 	-7.4925 	-0.14630 	0.10000 	0.27698 	38.18710;

# uniaxialMaterial 	Steel02 $matTag $Fy 	$E 			$b 			$R0 	$cR1 	$cR2 <$a1 $a2 $a3 $a4 $sigInit>
uniaxialMaterial 	Steel02 3 		75.0 	2805.72755 	0.01000 	18.0 	0.925 	0.15;
uniaxialMaterial 	Steel02 6 		75.0 	2805.72755 	0.01000 	18.0 	0.925 	0.15;

# Material for Column Shear
# uniaxialMaterial 	Elastic $matTag $E 				<$eta> <$Eneg>
uniaxialMaterial 	Elastic 7 		1000000000.0 	0.0;

# Material for flexural behavior
# limitCurve 	Rotation $curveTag 	$eleTag $dofl 	$dofv 	$iNodeTag 	$jNodeTag 	$fpc 	$fyt 	$Ag 	$rhot 	$thetay 	$VColOE 	$Kunload 		-$VyE
limitCurve 		Rotation 10000 		1001 	2 		3 		1 			4 			7.50 	75.0 	558.0 	0.004 	0.00 		340 		1002182.211525 	-VyE 	249;
limitCurve 		Rotation 10001 		1001 	1 		3 		1 			4 			7.50 	75.0 	558.0 	0.004 	0.00 		340 		1002182.211525 	-VyE 	249;

uniaxialMaterial LimitState 20000 100000.0 0.000998 1000000.0 0.010998 10000000.0 0.100998 -100000.0 -0.000998 -1000000.0 -0.010998 -10000000.0 -0.100998 0.810300 0.518800 0.0 0.0 0.688500 10000 3;
uniaxialMaterial LimitState 20001 100000.0 0.000998 1000000.0 0.010998 10000000.0 0.100998 -100000.0 -0.000998 -1000000.0 -0.010998 -10000000.0 -0.100998 0.810300 0.518800 0.0 0.0 0.688500 10001 3;

# limitCurve 	Rotation $curveTag 	$eleTag $dofl 	$dofv 	$iNodeTag 	$jNodeTag 	$fpc 	$fyt 	$Ag 	$rhot 	$thetay 	$VColOE 	$Kunload 		-$VyE
limitCurve 		Rotation 10002 		1004 	2 		3 		2 			7 			7.50 	75.0 	558.0 	0.004 	0.00 		340 		1002182.211525 	-VyE 	249;
limitCurve 		Rotation 10003 		1004 	1 		3 		2 			7 			7.50 	75.0 	558.0 	0.004	0.00 		340 		1002182.211525 	-VyE 	249;

uniaxialMaterial LimitState 20002 100000.0 0.000998 1000000.0 0.010998 10000000.0 0.100998 -100000.0 -0.000998 -1000000.0 -0.010998 -10000000.0 -0.100998 0.810300 0.518800 0.0 0.0 0.688500 10002 3;
uniaxialMaterial LimitState 20003 100000.0 0.000998 1000000.0 0.010998 10000000.0 0.100998 -100000.0 -0.000998 -1000000.0 -0.010998 -10000000.0 -0.100998 0.810300 0.518800 0.0 0.0 0.688500 10003 3;

############# Sections #############
# section 	Fiber $secTag 	<-GJ $GJ>
section 	Fiber 1 		-GJ 25370954.37936 {
	# patch quad $matTag 	$numSubdivIJ 	$numSubdivJK 	$yI 	$zI 	$yJ 	$zJ 	$yK 	$zK 	$yL 	$zL
	patch 	quad 1 			4 				8 				-12.9 	-10.7 	0.00 	-11.7 	0.00 	12.4 	-11.3 	11.1
	patch 	quad 1 			4 				8 				0.00 	-11.7 	12.9 	-10.7 	11.3 	11.1 	0.00 	12.4
	
	# layer straight $matTag 	$numFiber 	$areaFiber 	$yStart	$zStart $yEnd 	$zEnd
	layer 	straight 3 			2 			0.44 		0.00 	-8.70 	0.00 	9.10
	layer 	straight 3 			2 			1.00 		-9.70 	0.30 	9.70 	0.30
	layer 	straight 3 			2 			1.27 		-10.00 	-8.20 	10.00 	-8.20
	layer 	straight 3 			2 			1.27 		8.9 	8.70 	-8.9 	8.7
};

# section 	Fiber $secTag 	<-GJ $GJ>
section 	Fiber 2 		-GJ 25370954.37936 {
	# patch quad $matTag 	$numSubdivIJ 	$numSubdivJK 	$yI 	$zI 	$yJ 	$zJ 	$yK 	$zK 	$yL 	$zL
	patch 	quad 1 			4 				8 				-12.9 	-10.7 	0.00 	-11.7 	0.00 	12.4 	-11.3 	11.1
	patch 	quad 1 			4 				8 				0.00 	-11.7 	12.9 	-10.7 	11.3 	11.1 	0.00 	12.4
	
	# layer straight $matTag 	$numFiber 	$areaFiber 	$yStart	$zStart $yEnd 	$zEnd
	layer 	straight 3 			2 			0.44 		0.00 	-8.70 	0.00 	9.10
	layer 	straight 3 			2 			1.00 		-9.70 	0.30 	9.70 	0.30
	layer 	straight 3 			2 			1.27 		-10.00 	-8.20 	10.00 	-8.20
	layer 	straight 3 			2 			1.27 		8.9 	8.70 	-8.9 	8.7
};

# section Aggregator $secTag $matTag1 $dof1 $matTag2 $dof2 ....... <-section $sectionTag>
section Aggregator 3 7 Vy 7 Vz -section 1;
section Aggregator 4 7 Vy 7 Vz -section 2;

############# Geometric Transformations #############
geomTransf Linear 1 	0  1  0;

############# Elements #############
# element 	zeroLength $eleTag 	$iNode 	$jNode 	-mat $matTag1 	$matTag2 ... 						-dir $dir1 	$dir2 ... 	 	<-orient $x1 $x2 $x3 $yp1 $yp2 $yp3> 	<-doRayleigh $rFlag>
element 	zeroLength 1000 	1 		3 		-mat 7      	7      		7      7  20000  20001 	-dir 1  	2  3  4  5  6 	-orient 0.0 0.0 1.0 1.0 0.0 0.0;
element zeroLengthSection 1001 3 4 3 -orient 0.0 0.0 1.0 -1.0  0.0  0.0;

# element 	elasticBeamColumn $eleTag 	$iNode 	$jNode 	$A 		$E 			$G 			$J 			$Iy 		$Iz 		$transfTag <-mass $massDens> <-cMass>
element 	elasticBeamColumn 1002 		4 		5 		558.0 	493634.48	214623.687 	14784.120 	4994.272 	5481.24 	1; 		# FRAME ELEMENT FOR OFFSET w/ 100x E and G

# element 	zeroLength $eleTag 	$iNode 	$jNode 	-mat $matTag1 	$matTag2 ... 						-dir $dir1 	$dir2 ... 	 	<-orient $x1 $x2 $x3 $yp1 $yp2 $yp3> 	<-doRayleigh $rFlag>
element 	zeroLength 1003 	2 		6 		-mat 7      	7      		7      7  20002  20003 	-dir 1  	2  3  4  5  6 	-orient 0.0 0.0 1.0 1.0 0.0 0.0;
element zeroLengthSection 1004 6 7 4 -orient 0.0 0.0 -1.0 -1.0  0.0  0.0;

# element 	elasticBeamColumn $eleTag 	$iNode 	$jNode 	$A 		$E 			$G 			$J 			$Iy 		$Iz 		$transfTag <-mass $massDens> <-cMass>
element 	elasticBeamColumn 1005 		7 		8 		558.000 493634.48 	214623.687 	14784.120 	4994.272 	5481.24		1; 		# FRAME ELEMENT FOR OFFSET w/ 100x E and G
element 	elasticBeamColumn 1006 		5 		8 		558.000 4936.3448 	2146.23687 	14784.120 	4994.272 	5481.24 	1; 		# MAIN FRAME ELEMENT

############# Recorders #############
recorder Node 		-file ./out/1000C_2497kGravity_Cyclic/node_disp.out 		-node 1 2 3 4 5 6 7 8 	-dof 1 disp;
recorder Node 		-file ./out/1000C_2497kGravity_Cyclic/col_force.out 		-node 1 				-dof 1 reaction;
recorder Node 		-file ./out/1000C_2497kGravity_Cyclic/rotation_top.out 		-node 1 2 3 4 5 6 7 8 	-dof 5 disp;
recorder Element 	-file ./out/1000C_2497kGravity_Cyclic/fail.out -time 		-ele 1000 material 6 failFlag;
recorder Element 	-file ./out/1000C_2497kGravity_Cyclic/force_fiber.out 		-ele 1004 force;
recorder Element 	-file ./out/1000C_2497kGravity_Cyclic/deformation_fiber.out -ele 1001 deformation;
recorder Element 	-file ./out/1000C_2497kGravity_Cyclic/rotation_j.out 		-ele 1003 material 6 rotationAndLimits;
record;

####################################################################################
## Analysis
####################################################################################

####################################################################################
## Analysis 1 - GRAVITY
####################################################################################
timeSeries Linear 1;
pattern Plain 1 1 {
		# load 2 0.0    0.0    	0.0    0.0    0.0    0.0
		# load 2 0.0    0.0    -113.0    0.0    0.0    0.0
		# load 2 0.0    0.0    -155.4    0.0    0.0    0.0 	
		# load 2 0.0    0.0    -466.2    0.0    0.0    0.0
		# load 2 0.0    0.0    -621.6    0.0    0.0    0.0
		
		# load 2 0.0    0.0    -416.0    0.0    0.0    0.0
		# load 2 0.0    0.0    -1249.0    0.0    0.0    0.0
	load 2 0.0    0.0    -2497.0    0.0    0.0    0.0
};

constraints Plain;
numberer RCM;
system SparseGEN;
test NormDispIncr 1.0e-8 600;
algorithm Newton;
integrator LoadControl 0.1;
analysis Static;
analyze 10;
loadConst -time 0.0;
puts "\n> gravity load applied";
wipeAnalysis;

####################################################################################
## Analysis 2 - (CYCLIC) PUSHOVER
####################################################################################
timeSeries Linear 2;
pattern Plain 2 2 {
	load 2 1.00000 0.00000 0.00000 0.00000 0.00000 0.00000
};

set ctrlNode 2;
set ctrlDOF 1;
set targIncr 1e-4;
# set targTol 0.00000001;
set targTol 1.0e-8;
set maxIter 10000;

# foreach targDisp {0           0           0           0           0        0.02        0.05         0.1        0.14        0.15        0.12         0.1        0.07        0.02           0       -0.05       -0.08       -0.15       -0.12       -0.12       -0.12       -0.12       -0.12        -0.1       -0.07       -0.05       -0.03           0        0.02        0.04        0.08         0.1        0.12        0.15        0.12        0.09        0.07        0.04        0.02           0       -0.02       -0.05       -0.08        -0.1       -0.13       -0.15       -0.12       -0.09       -0.07       -0.05       -0.03           0        0.02        0.04        0.07        0.09        0.11        0.14        0.12         0.1        0.07        0.04        0.02       -0.01       -0.03       -0.05       -0.08       -0.11       -0.13       -0.15       -0.12        -0.1       -0.07       -0.05       -0.03           0        0.02        0.05        0.07        0.14         0.2        0.25         0.3        0.24         0.2        0.15        0.09        0.04       -0.01       -0.05        -0.1       -0.15       -0.21       -0.25        -0.3       -0.25        -0.2       -0.15        -0.1       -0.06           0        0.04        0.09        0.15        0.19        0.24        0.29        0.24        0.19        0.15        0.09        0.04       -0.01       -0.05        -0.1       -0.15       -0.21       -0.25       -0.31       -0.25       -0.21       -0.15        -0.1       -0.05       -0.01        0.05         0.1        0.14         0.2        0.24         0.3        0.25         0.2        0.14         0.1        0.05       -0.01       -0.06        -0.1       -0.15       -0.21       -0.25        -0.3       -0.25       -0.21       -0.15        -0.1       -0.05           0        0.09         0.2         0.3        0.39        0.49        0.59        0.39        0.29        0.19         0.1       -0.01        -0.1       -0.21       -0.31        -0.4        -0.5        -0.6       -0.51       -0.34        -0.3       -0.21        -0.1       -0.01         0.1        0.19         0.3        0.39        0.49        0.59        0.39        0.29        0.19        0.09       -0.01        -0.1        -0.2        -0.3        -0.4       -0.51        -0.6       -0.51       -0.41        -0.3       -0.21        -0.1           0        0.13        0.19         0.3        0.39        0.49        0.59        0.49        0.39        0.29         0.2        0.09       -0.01        -0.1       -0.21        -0.3        -0.4        -0.5        -0.6        -0.5        -0.4        -0.3       -0.21        -0.1         0.2         0.4        0.59        0.84        0.95        1.08        1.04        0.94        0.74         0.6        0.39        0.29        0.27        0.09        -0.1       -0.31        -0.5       -0.69       -0.83       -0.94       -1.07        -1.2       -1.17       -1.08       -0.87       -0.69       -0.48       -0.28       -0.09        0.11        0.31         0.5        0.71        0.93        1.17        0.98         0.7        0.63        0.43        0.24        0.07           0       -0.12       -0.31       -0.53       -0.71       -0.91       -1.09        -1.2       -1.09       -0.88       -0.68       -0.48        -0.3       -0.08        0.11        0.31         0.5         0.7         0.9        1.09        1.18        1.05        0.88        0.68        0.49        0.28        0.13       -0.07       -0.27       -0.41       -0.52       -0.71        -0.9       -1.09        -1.2       -1.07       -0.88       -0.71        -0.5       -0.31       -0.09        0.11        0.31        0.51         0.7        0.94        1.12        1.24        1.38        1.54        1.65        1.75         1.8        1.74        1.51        1.32        1.13        0.94         0.8        0.63        0.43        0.22        0.01       -0.16       -0.36       -0.52       -0.71       -0.92       -1.12       -1.33       -1.51       -1.69       -1.81       -1.74       -1.38       -1.24       -1.04       -0.84       -0.64        -0.4        -0.2        0.02         0.2        0.41        0.59        0.79        0.99        1.17        1.35        1.58        1.69         1.8        1.69        1.36        0.98        0.34        0.18       -0.02       -0.21       -0.41       -0.61       -0.81       -1.02       -1.18       -1.42       -1.62        -1.8       -1.68       -1.49       -1.28       -1.02       -0.79       -0.58       -0.36        -0.1        0.14        0.46        0.77        0.96        1.16        1.35        1.59        1.79        1.69        1.53        1.32        1.09        0.94        0.76        0.45        0.08       -0.22       -0.49       -0.71       -0.91       -1.13       -1.31       -1.52       -1.69       -1.81       -1.68       -1.45       -1.22       -1.04       -0.84       -0.49       -0.29       -0.07        0.11        0.31         0.5         0.7        0.87        1.02        1.24        1.46        1.62        1.81        1.91        2.06        2.23         2.4        2.35        2.13        1.93        1.74        1.53        1.33        1.14        0.93        0.73        0.51        0.33        0.13       -0.07       -0.27       -0.45       -0.65       -0.86       -1.07       -1.25       -1.44       -1.68       -1.88       -2.09       -2.28       -2.41       -2.31        -2.1        -1.9       -1.69        -1.5       -1.28       -1.08       -0.87       -0.69        -0.5        -0.3       -0.06        0.12        0.31        0.51        0.69        0.89        1.11        1.31         1.5        1.66        1.91        2.12        2.28        2.39        2.24        2.07        1.89         1.7         1.5        1.34        1.13        0.92        0.72        0.51        0.32        0.14       -0.04       -0.22       -0.42       -0.61       -0.77       -0.91       -1.06        -1.2       -1.35        -1.5       -1.63       -1.76       -1.91       -2.13       -2.28       -2.41       -2.18          -2        -1.8       -1.63        -1.4       -1.21       -1.02       -0.83       -0.61       -0.37       -0.11        0.11        0.37        0.61         0.9        1.13        1.37        1.61        1.85        2.09         2.3        2.39        2.31        2.16        1.96        1.76        1.57        1.39         1.2        0.94         0.7         0.5        0.27        0.03       -0.18       -0.42       -0.64       -0.87       -1.06       -1.27       -1.46       -1.66       -1.86       -2.07       -2.27       -2.41       -2.29       -2.18       -2.03       -1.88       -1.72       -1.58       -1.41       -1.13       -1.12       -1.06       -0.48       -0.18        0.02        0.25         0.4         0.6        0.77        0.94        1.13        1.32        1.57        1.84         2.1        2.43        2.58        2.72        2.84        2.98           3        2.89        2.76         2.6        2.23         2.1        1.86        1.65        1.27        1.09        0.88         0.6        0.33       -0.02       -0.45       -0.86       -1.28       -1.61       -1.94        -2.2       -2.47       -2.68        -2.9       -3.01       -2.98        -2.8       -2.33       -2.17       -1.96       -1.49       -1.16       -0.83       -0.49        -0.2        0.11        0.48        0.83        1.13         1.4        1.69        1.97         2.2        2.39        2.57        2.66        2.98         2.9        2.72        2.49        2.12        1.67        1.25        0.76        0.28       -0.21       -0.79       -1.24       -1.64       -2.13       -2.94       -3.01       -2.52        -2.2       -1.81        -1.5       -1.19       -0.69       -0.28        0.18        0.67           1        1.43        1.86        2.18        2.45        2.66           3        3.01        2.94        2.72        2.46        1.83        1.23         0.5       -0.02       -0.66       -1.34       -1.61       -1.82       -1.83       -1.88       -1.89       -1.97       -2.09       -2.18       -2.28        -2.4        -2.5        -2.6       -2.71       -2.83       -2.98       -3.02       -2.95       -2.87       -2.82       -2.23       -2.01       -1.21       -0.72           0        0.76        1.16        1.33        1.56        1.75        1.94        2.39        2.69        2.87        3.02        3.14        3.18         3.2        3.27        3.33        3.34        3.34        3.35        3.37        3.42        3.43        3.43        3.43        3.43        3.43        3.43        3.43        3.43        3.43        3.43        3.53        3.58        3.64         3.6         3.5        3.39        3.13        2.79        2.49        2.15        2.02        1.62        1.28        0.93        0.32       -0.08       -0.25       -0.59       -0.91       -1.29       -1.66        -0.3} 
# foreach targDisp {0.9525 	 	-0.9525 	 	0.9525 	 	-0.9525 	 	1.4325 	 	-1.4325 	 	1.4325 	 	-1.4325 	 	2.1525 	 	-2.1525 	 	2.1525 	 	-2.1525 	 	3.0075 	 	-3.0075 	 	3.0075 	 	-3.0075 	 	4.215 	 	-4.215 	 	4.215 	 	-4.215 	 	5.9025 	 	-5.9025 	 	5.9025 	 	-5.9025 	 	6.6075 	 	-6.6075 	 	6.6075 	 	-6.6075 	 	7.5975 	 	-7.5975 	 	7.5975 	 	-7.5975 	 	8.7375 	 	-8.7375 	 	8.7375 	 	-8.7375 	 	10.05 	 	-10.05 	 	10.05 	 	-10.05} {
foreach targDisp {0            6.0} {

	constraints Plain;
	numberer RCM;
	system Umfpack;
	set ctrlDisp [nodeDisp $ctrlNode $ctrlDOF];
	set travel 0.0;
	set relDisp [expr $targDisp - $ctrlDisp];
	if {$relDisp > 0} {
		set sgn 1.0;
	} else {
		set sgn -1.0;
	};
	set incr [expr $sgn*$targIncr];
	puts "> excursion: $targDisp | increment: $incr";
	while {[expr abs($travel)] < [expr abs($relDisp)]} {
		test NormDispIncr $targTol $maxIter;
		algorithm Newton;
		integrator DisplacementControl $ctrlNode $ctrlDOF $incr;
		analysis Static;
		set ok [analyze 1];
		if {$ok == 0} {
			set travel [expr $travel + $incr];
		} elseif {$ok != 0} {
			set prntDisp [expr int([nodeDisp $ctrlNode $ctrlDOF]*100.0)/100.0];
			puts "	> at $prntDisp";
			set tempIncr $incr;
			set denom 10.0;
			set counter 0;
			set tempTol $targTol;
			while {$ok != 0} {
				incr counter;
				if {$counter == 1} {
					algorithm KrylovNewton -maxDim 4;
				} elseif {$counter == 10} {
					exit;
				};
				set tempIncr [expr $tempIncr/$denom];
				puts "		> trying increment: $tempIncr (counter = $counter)\n";
				integrator DisplacementControl $ctrlNode $ctrlDOF $tempIncr;
				set ok [analyze 1];
			};
			set travel [expr $travel + $tempIncr];
		};
	};
};

puts "Analysis of Model 1000C_2497kGravity_Cyclic Complete"
wipe;
#exit;

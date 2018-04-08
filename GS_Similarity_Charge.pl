##!/usr/bin/env perl





use strict;
use warnings;
# install: 
#          sudo cpan Parallel::ForkManager
use Parallel::ForkManager;
use Benchmark;




#
# Global variable
my $numb_atoms;
#
my %Atomic_number = ( '89'  => 'Ac', '13'  => 'Al', '95'  => 'Am', '51'  => 'Sb',	
	                  '18'  => 'Ar', '33'  => 'As', '85'  => 'At', '16'  => 'S',  
					  '56'  => 'Ba', '4'   => 'Be', '97'  => 'Bk', '83'  => 'Bi',	
                      '107' => 'Bh', '5'   => 'B', 	'35'  => 'Br', '48'  => 'Cd',	
	                  '20'  => 'Ca', '98'  => 'Cf',	'6'   => 'C',  '58'  => 'Ce',	
	                  '55'  => 'Cs', '17'  => 'Cl',	'27'  => 'Co', '29'  => 'Cu',	
	                  '24'  => 'Cr', '96'  => 'Cm', '110' => 'Ds', '66'  => 'Dy',	
	                  '105' => 'Db', '99'  => 'Es', '68'  => 'Er', '21'  => 'Sc',	
	                  '50'  => 'Sn', '38'  => 'Sr', '63'  => 'Eu', '100' => 'Fm',	
	                  '9'   => 'F',  '15'  => 'P',  '87'  => 'Fr', '64'  => 'Gd',	
	                  '31'  => 'Ga', '32'  => 'Ge', '72'  => 'Hf', '108' => 'Hs',	
                      '2'   => 'He', '1'   => 'H',  '26'  => 'Fe', '67'  => 'Ho',	
					  '49'  => 'In', '53'  => 'I',  '77'  => 'Ir', '70'  => 'Yb',
					  '39'  => 'Y',  '36'  => 'Kr', '57'  => 'La', '103' => 'Lr',	
					  '3'   => 'Li', '71'  => 'Lu', '12'  => 'Mg', '25'  => 'Mn',	
                      '109' => 'Mt', '101' => 'Md', '80'  => 'Hg', '42'  => 'Mo',	
					  '60'  => 'Nd', '10'  => 'Ne', '93'  => 'Np', '41'  => 'Nb',	
					  '28'  => 'Ni', '7'   => 'N',  '102' => 'No', '79'  => 'Au',	
					  '76'  => 'Os', '8'   => 'O', 	'46'  => 'Pd', '47'  => 'Ag',	
					  '78'  => 'Pt', '82'  => 'Pb',	'94'  => 'Pu', '84'  => 'Po',	
					  '19'  => 'K',  '59'  => 'Pr', '61'  => 'Pm', '91'  => 'Pa',	
					  '88'  => 'Ra', '86'  => 'Rn', '75'  => 'Re', '45'  => 'Rh',	
					  '37'  => 'Rb', '44'  => 'Ru', '104' => 'Rf', '62'  => 'Sm',
					  '106' => 'Sg', '34'  => 'Se', '14'  => 'Si', '11'  => 'Na',
					  '81'  => 'Tl', '73'  => 'Ta', '43'  => 'Tc', '52'  => 'Te',	
					  '65'  => 'Tb', '22'  => 'Ti', '90'  => 'Th', '69'  => 'Tm',	
					  '112' => 'Uub','116' => 'Uuh','111' => 'Uuu','118' => 'Uuo',	
					  '115' => 'Uup','114' => 'Uuq','117' => 'Uus','113' => 'Uut',
					  '92'  => 'U',  '23'  => 'V',  '74'  => 'W',  '54'  => 'Xe',
                      '30'  => 'Zn', '40'  => 'Zr' );
					  
###################################
# read files
sub read_file {
	# filename
	my ($input_file) = @_;
	my @lines = ();
	my @array = ();
	# open file
	open(FILE, "<", $input_file ) || die "Can't open $input_file: $!";
	@lines = <FILE>;
	close (FILE);
	# loop
	foreach $a (@lines){
		chomp ($a);
		$array[++$#array] = $a;
	}
	# return array	
	return @array;
}
###################################
# Cartesian Coords
sub cartesian_coords {
	# filename
	my ($input_file) = @_;
	my @matrix       = read_file($input_file);
	my @total_coords     = ();
	my @mulliken_charges = ();
	my @APT_charges      = ();
	# variable
	my $num_atoms_global;
	my $global_energy;
	my $dipole_moment;
	# coodenadas
	my @coords;
	my @columns_1N = ();
	my @columns_2N = ();
	my @columns_3N = ();
	my @columns_4N = ();
	my @columns_5N = ();
	my @columns_6N = ();
	#
	my $count_lines = 0;
	foreach my $a_1 (@matrix){
		# SCF Done:  E(RPBE1PBE) =  -56.7829127857     A.U. after   40 cycles
		if ( ($a_1=~/SCF/gi ) && ($a_1=~/Done/gi ) && ($a_1=~/after/gi ) ){
			my @array_tabs = ();
			@array_tabs    = split (/ /,$a_1);
			push (@columns_1N,$array_tabs[7]);
		}
		# Standard orientation:
		if ( ($a_1=~/Standard/gi ) && ($a_1=~/orientation/gi ) && ($a_1=~/:/gi ) ){
			push (@columns_2N,$count_lines);
		}
		# Rotational constants (GHZ):
		if ( ($a_1=~/Rotational/gi ) && ($a_1=~/constants/gi ) && ($a_1=~/GHZ/gi ) ){
			push (@columns_3N,$count_lines);
		}
		# Dipole moment (field-independent basis, Debye):
		# X= 0.0001 Y= -0.0003 Z= 0.8792 Tot= 0.8792
		if ( ($a_1=~/Dipole/gi ) && ($a_1=~/moment/gi ) && ($a_1=~/basis/gi )  && ($a_1=~/Debye/gi ) ){
			my $sumTmp = $count_lines + 1;
			push (@columns_4N,$sumTmp);
		}
		# Mulliken charges:
		if ( ($a_1=~/Mulliken/gi ) && ($a_1=~/charges:/gi ) ){
			my $sumTmp = $count_lines + 1;
			push (@columns_5N,$sumTmp);
		}		
		# APT charges:
		if ( ($a_1=~/APT/gi ) && ($a_1=~/charges:/gi ) ){
			my $sumTmp = $count_lines + 1;
			push (@columns_6N,$sumTmp);
		}
		$count_lines++;
	}
	if ( scalar (@columns_1N) > 0 ){
		for (my $i=0; $i < scalar (@columns_1N); $i++){
			#
			my $start         = $columns_2N[$i] + 5;
			my $end           = $columns_3N[$i] - 2;
			$num_atoms_global = $end - $start + 1;
			$numb_atoms       = $num_atoms_global;
			#
			$global_energy = $columns_1N[$i];
			@coords = ();
			foreach my $j (@matrix[$start..$end]){
				push (@coords,$j);		
			}
		}
		foreach my $i (@coords){
			my @tmp = ();
			@tmp    = split (/\s+/,$i);
			push (@total_coords,"$tmp[2]\t$tmp[4]\t$tmp[5]\t$tmp[6]");
		}				
		my $line = (scalar(@columns_4N)) - 2;
		my $data = $columns_4N[$line];
		my @aux  = ();
		@aux     = split (/\s+/,$matrix[$data]);
		$dipole_moment = "$aux[2]\t$aux[4]\t$aux[6]";
		#
		my $leng1   = scalar (@columns_5N);
		my $da1     = $leng1 - 1;
		my $inicio1 = $columns_5N[$da1] + 1 ;
		my $final1  = $inicio1 + $numb_atoms - 1 ;
		foreach my $j (@matrix[$inicio1..$final1]){
			my @MyLa = split (/\s+/,$j);
			push (@mulliken_charges,$MyLa[3]);
		}
		#
		my $leng2   = scalar (@columns_6N);
		my $da2     = $leng2 - 1;
		my $inicio2 = $columns_6N[$da2] + 1 ;
		my $final2  = $inicio2 + $numb_atoms - 1 ;
		foreach my $j (@matrix[$inicio2..$final2]){
			my @MyLa = split (/\s+/,$j);
			push (@APT_charges,$MyLa[3]);
		}
		# array content array
		my @array_set  = ([@total_coords],
						  [$num_atoms_global],
                          [$global_energy],
                          [$dipole_moment],
						  [@mulliken_charges],
						  [@APT_charges]);

		return @array_set; 
	} else {
		print "No Present SCF\n";
	}
}
###################################
# promedio
sub promedio {
	my ($num,$data) = @_;
	# write file
	my $sum = 0;
	for ( my $i = 0 ; $i < $num ; $i = $i + 1 ){
		$sum+= @$data[$i];
	}
	my $div = $sum / $num;
	return $div; 
}
###################################
# delete repeat data
sub uniq {
	my %seen;
	grep !$seen{$_}++, @_;
}
###################################
# index duplicate data
sub index_elements {
	my ($duplicate_name,$files_name) = @_;
	# reference arrays	
	my @array_1     = @{$duplicate_name}; 
	my @array_2     = @{$files_name};
	my @array_index = ();
	#
	my @filtered = uniq(@array_1);
	foreach my $u (@filtered){
		my @del_indexes = reverse( grep { $array_2[$_] eq "$u" } 0..$#array_2);
		foreach my $k (@del_indexes) {
			push (@array_index,$k);
		}
	}
	return @array_index;
}

###################################
# MAIN
#
my $tiempo_inicial  = new Benchmark; #funcion para el tiempo de ejecucion del programa

my ($threshold_duplicate) = @ARGV;
if (not defined $threshold_duplicate) {
	die "\nGrigoryan Springborg Similarity Charge must be run with:\n\nUsage:\n\tGS_Similarity_Charge.pl [threshold duplicate]\n";
	exit(1);  
}

#my @files_g09   = glob "./Outs/*.log";
my @files_g09   = glob "./Outs/*.out";

my $ncpus = 100;
my $pm    = new Parallel::ForkManager($ncpus);
my $iteration = 0;
# Array elements and coords
my @atoms_coords  = ();
my @array_keys    = ();
my @files_sn_exte = ();
#
my %Info_Charge   = ();
my %Info_Coords   = ();
my %Info_Energy   = ();
my %Info_Name     = ();
#
for ( my $x = 0 ; $x < scalar (@files_g09) ; $x = $x + 1 ){
	my $file_name          = $files_g09[$x];
	(my $without_extension = $file_name) =~ s/\.[^.]+$//;
	push (@files_sn_exte,$without_extension);
	#
	my @matrix           = cartesian_coords($file_name);
	my $element_coords;
	for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
		my @array_tabs  = split (/\s+/,$matrix[0][$i]);
		my $radii_val;
		my $other_element = 0;
		if ( exists $Atomic_number{$array_tabs[0]} ) {
			# Exists
			$radii_val = $Atomic_number{$array_tabs[0]};
		} else {
			# Not exists
			$radii_val = $other_element ;
		}
		$element_coords.= "$radii_val  $array_tabs[1]  $array_tabs[2]  $array_tabs[3]\n";
	}
	push (@atoms_coords,$element_coords);
	#
	my @Charges_Mol = ();
	for ( my $i = 0 ; $i < $numb_atoms ; $i = $i + 1 ){
		for ( my $j = 0 ; $j < $numb_atoms ; $j = $j + 1 ){
			if ( $i < $j ){
				my $mult_1 = $matrix[4][$i] * $matrix[4][$j];
				my $mult_2 = $matrix[4][$i] * $matrix[4][$j];
				push (@Charges_Mol,$mult_1);
			}
		}
	}
	#
	my $id            = sprintf("%06d",$x);
	#
	$Info_Energy{$id} = $matrix[2][0];
	$Info_Coords{$id} = $element_coords;
	$Info_Name{$id}   = $without_extension;
	$Info_Charge{$id} = \@Charges_Mol;
	push (@array_keys,$id);
}
#
#
my $file_tmp = "Dupli.tmp";
open (FILE, ">$file_tmp") or die "Unable to open file: $file_tmp";
my $file_log = "Info_Duplicates.txt";
open (LOGDUPLI, ">$file_log") or die "Unable to open file: $file_log"; 
print LOGDUPLI "\n# # # SUMMARY SIMILAR STRUCTURES # # #\n\n";
#
for ( my $x = 0 ; $x < scalar (@files_sn_exte) ; $x = $x + 1 ){
	$pm->start($iteration) and next;
	# All children process havee their own random.			
	srand();
	for ( my $y = 0 ; $y < scalar (@files_sn_exte) ; $y = $y + 1 ){
		if ( $x < $y ){
			my @Charges_Mol1 = @{$Info_Charge{$array_keys[$x]}};
			my @Charges_Mol2 = @{$Info_Charge{$array_keys[$y]}};
			#
			my $InterDist_1 = (2/($numb_atoms*($numb_atoms-1)));
			my $InterDist_2 = (($numb_atoms*($numb_atoms-1))/2);  
			#
			my @mol_alpha = ();
			my @mol_beta  = ();
			my @idx_1 = sort { $Charges_Mol1[$a] <=> $Charges_Mol1[$b] } 0 .. $#Charges_Mol1;
			my @idx_2 = sort { $Charges_Mol2[$a] <=> $Charges_Mol2[$b] } 0 .. $#Charges_Mol2;
			@mol_alpha = @Charges_Mol1[@idx_1];
			@mol_beta  = @Charges_Mol2[@idx_2];
			#
			my $num_1 = scalar (@mol_alpha);
			my $num_2 = scalar (@mol_beta);
			my $dim_alpha =  promedio ($num_1,\@mol_alpha);
			my $dim_beta  =  promedio ($num_2,\@mol_beta);
			#
			my $sumX;
			my $sumY;
			# Sin normalizar
			for ( my $i = 0 ; $i < $InterDist_2 ; $i = $i + 1 ){
				my $mult = ( $mol_alpha[$i] - $mol_beta[$i] )**2;
				$sumX+=$mult;
			}
			my $Springborg_1 = sqrt( $InterDist_1 * $sumX );
			# Normalizado
			for ( my $i = 0 ; $i < $InterDist_2 ; $i = $i + 1 ){
				my $mult = ( ($mol_alpha[$i]/$dim_alpha) - ($mol_beta[$i]/$dim_beta) )**2;
				$sumY+=$mult;
			}
			my $Springborg_2 = sqrt( $InterDist_1 * $sumY );
			#
			if ( $Springborg_2 < $threshold_duplicate ) {			
				my $number      = sprintf '%.6f', $Springborg_2;
				print FILE "$array_keys[$y]\n";
				print FILE "Value = $number\n";
				print LOGDUPLI "# $files_sn_exte[$x] ~= $files_sn_exte[$y]\n";
				print LOGDUPLI "# Value = $number\n";
				print LOGDUPLI "------------------------\n";
			}
		}
		$iteration++;
	}
	$pm->finish;	
}
close (FILE);
# Paralel
$pm->wait_all_children;
# # #
my @data_tmp = read_file ($file_tmp);
my @duplicates_name = ();
my @Value_simi      = ();
foreach my $info (@data_tmp) {
	if ( ($info =~ m/Value/) ) {
		my @array_tabs = ();
		@array_tabs    = split ('\s+',$info);
		push (@Value_simi,$array_tabs[2]);
	} else {
		push (@duplicates_name,$info);
	}
}
# Delete similar structures
my @index_files = index_elements (\@duplicates_name,\@array_keys);

my $file_xyz = "02Duplicates_coords.xyz";
open (DUPLIXYZ, ">$file_xyz") or die "Unable to open file: $file_log"; 
my $count_sim_struc = 0;
foreach my $id_index (@index_files) {
	my $coords_dup = $Info_Coords{$array_keys[$id_index]};
	print DUPLIXYZ "$numb_atoms\n";
	print DUPLIXYZ "Duplicate $files_sn_exte[$id_index]\n";
	print DUPLIXYZ "$coords_dup";
	$count_sim_struc++;
}
print LOGDUPLI "\nNumber of Similar Structures = $count_sim_struc\n";
close (LOGDUPLI);  
close (DUPLIXYZ);
#
# Delete similar structures
for my $k (@index_files) {
	delete $Info_Coords{$array_keys[$k]};
	delete $Info_Energy{$array_keys[$k]};
}
#
my @energy_isomer = ();
my @coords_isomer = ();
my @name_isomer   = ();
foreach my $key (sort(keys %Info_Coords)) {
	push (@name_isomer  ,$Info_Name{$key});
	push (@energy_isomer,$Info_Energy{$key});
	push (@coords_isomer,$Info_Coords{$key});
}

my @mol_energy = ();
my @mol_coords = ();
my @idx        = sort { $energy_isomer[$a] <=> $energy_isomer[$b] } 0 .. $#energy_isomer;
@mol_energy    = @energy_isomer[@idx];
@mol_coords    = @coords_isomer[@idx];
#
my $file_x = "01Clean_Duplicates_coords.xyz";
open (CLEANDUPLIXYZ, ">$file_x") or die "Unable to open XYZ file: $file_x"; 
for ( my $x = 0 ; $x < scalar (@mol_energy) ; $x = $x + 1 ){
	# resta
	my $resta = abs($mol_energy[0]) - abs($mol_energy[$x]);
	# 1 Hartree = 27,2114 ev
	# 1 Hartree = 627,509 Kcal/mol
	my $eV      = sprintf("%06f",(27.2114 * $resta ));
	my $Kcalmol = sprintf("%06f",(627.509 * $resta ));
	my $Hartree = sprintf("%06f",$mol_energy[$x]);
	print CLEANDUPLIXYZ "$numb_atoms\n";
	print CLEANDUPLIXYZ "$Kcalmol Kcal/mol $eV eV $Hartree H $name_isomer[$x]\n";
	print CLEANDUPLIXYZ "$mol_coords[$x]";
}
close (CLEANDUPLIXYZ);
#
my $tiempo_final  = new Benchmark;
my $tiempo_total  = timediff($tiempo_final, $tiempo_inicial);
print "\n\tExecution Time: ",timestr($tiempo_total),"\n";
print "\n";

#
unlink ($file_tmp);
exit 0;

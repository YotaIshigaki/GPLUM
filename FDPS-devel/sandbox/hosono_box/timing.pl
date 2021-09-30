use vars qw(@DATA);
use strict;

sub ReadFile{
	my $file = $_[0];
	local $/ = "";
	open(DAT, $file);
	my @data = <DAT>;
	close DAT;
	return @data;
}

sub ParseData{
	my $data = $_[0];
	#alloc local vars
	my %cnt;
	my %max_time;
	my %result;
	#
	foreach my $line (split/\n/, $data){
		if($line =~ /^([\w\s]+):/){
			my $key = $1;
			$cnt{"$key"} += 1;
			if($line =~ /max_time= ([\d\.\-\+e]+)/){
				$max_time{$key . $cnt{$key}} = $1;
			}
		}
		if($line =~ /total_speed=([\d\.\-\+e]+)(.+)wtime_tot=([\d\.\-\+e]+)/){
			$result{"Tflops"} = $1;
			$result{"tottime"} = $3;
		}
	}
	############
	$result{"exPtcl"} += $max_time{"exchangeParticle1"};
	$result{"decDom"} += $max_time{"decomposeDomainMultiStep1"};
	my $loop = $cnt{"setParticleLocalTree"} - 2;
	for(my $i = 1 ; $i <= $loop ; $i += 1){
		$result{"DensDrvtExLET"} += ($max_time{"exchangeLocalEssentialTree$i"});
		$result{"DensDrvtClcFo"} += ($max_time{"calcForce$i"});
	}	

	$result{"HydrExLET"} += ($max_time{"exchangeLocalEssentialTree" . ($loop+1) });
	$result{"HydrClcFo"} += ($max_time{"calcForce" . ($loop+1)});

	$result{"GravExLET"} += ($max_time{"exchangeLocalEssentialTree" . ($loop+2) });
	$result{"GravClcFo"} += ($max_time{"calcForce" . ($loop+2)});


	foreach my $key(keys %result){
		print "$key = " . $result{$key} . "\n"
	}
	#foreach my $key(keys %max_time){
	#	print "$key = " . $max_time{$key} . "\n"
	#}

}

@DATA = ReadFile("time_0000.dat");
ParseData($DATA[$#DATA-1]);





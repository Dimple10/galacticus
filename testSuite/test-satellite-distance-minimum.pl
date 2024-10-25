#!/usr/bin/env perl
use strict;
use warnings;
use lib $ENV{'GALACTICUS_EXEC_PATH'         }."/perl";
use lib $ENV{'GALACTICUS_ANALYSIS_PERL_PATH'}."/perl";
use PDL;
use PDL::NiceSlice;
use PDL::IO::HDF5;
use File::Slurp qw(slurp);
use System::Redirect;

# Check calculations of minimum satellite distance from host center.
# Andrew Benson (28-November-2022)

# Run the model.
&System::Redirect::tofile("cd ..; ./Galacticus.exe testSuite/parameters/satelliteDistanceMinimum.xml","outputs/satelliteDistanceMinimum.log");
unless ( $? == 0 ) {
    print "FAILED:  model run:\n";
    system("cat outputs/satelliteDistanceMinimum.log");
} else {
    print "SUCCESS: model run\n";
}

# Read the model data and check for consistency.
my $model   = new PDL::IO::HDF5("outputs/satelliteDistanceMinimum.hdf5");
my $outputs = $model->group('Outputs');
foreach my $outputName ( $outputs->groups() ) {
    my $output            = $outputs ->group  ($outputName               )       ;
    my $nodeData          = $output  ->group  ('nodeData'                )       ;
    my $isIsolated        = $nodeData->dataset('nodeIsIsolated'          )->get();
    my $nodeIndex         = $nodeData->dataset('nodeIndex'               )->get();
    my $distanceMinimum   = $nodeData->dataset('satelliteDistanceMinimum')->get();
    my $positionX         = $nodeData->dataset('satellitePositionX'      )->get();
    my $positionY         = $nodeData->dataset('satellitePositionY'      )->get();
    my $positionZ         = $nodeData->dataset('satellitePositionZ'      )->get();
    my $distance          =
	sqrt
	(
	 +$positionX**2
	 +$positionY**2
	 +$positionZ**2
	);
    my $satellites        =
	which
	(
	 $isIsolated == 0
	);
    my $tolerance = 1.0e-6;
    my $status    = all($distanceMinimum->($satellites) <= $distance->($satellites)*(1.0+$tolerance)) ? "SUCCESS" : "FAILED";
    my $offset    = $distanceMinimum->($satellites)/$distance->($satellites)-1.0;
    print "Maximum fractional offset: ".$offset->maximum()."\n";
    print $status.": minimum distance from host center\n";
}

exit 0;

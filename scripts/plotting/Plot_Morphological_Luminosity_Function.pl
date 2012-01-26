#!/usr/bin/env perl
my $galacticusPath;
if ( exists($ENV{"GALACTICUS_ROOT_V091"}) ) {
 $galacticusPath = $ENV{"GALACTICUS_ROOT_V091"};
 $galacticusPath .= "/" unless ( $galacticusPath =~ m/\/$/ );
} else {
 $galacticusPath = "./";
}
unshift(@INC,$galacticusPath."perl"); 
use PDL;
use PDL::NiceSlice;
use XML::Simple;
use Math::SigFigs;
use Data::Dumper;
require Galacticus::HDF5;
require Galacticus::Magnitudes;
require Stats::Histograms;
require GnuPlot::PrettyPlots;
require GnuPlot::LaTeX;
require XMP::MetaData;

# Get name of input and output files.
if ( $#ARGV != 1 && $#ARGV != 2 ) {die("Plot_Morphological_Luminosity_Function.pl <galacticusFile> <outputDir/File> [<showFit>]")};
$self           = $0;
$galacticusFile = $ARGV[0];
$outputTo       = $ARGV[1];
if ( $#ARGV == 2 ) {
    $showFit    = $ARGV[2];
    if ( lc($showFit) eq "showfit"   ) {$showFit = 1};
    if ( lc($showFit) eq "noshowfit" ) {$showFit = 0};
} else {
    $showFit = 0;
}

# Check if output location is file or directory.
my $outputDir;
if ( $outputTo =~ m/\.pdf$/ ) {
    $outputFile = $outputTo;
    $outputDir = ".";
} else {
    system("mkdir -p $outputTo");
    $outputFile = $outputTo."/Morphological_Luminosity_Function.pdf";
    $outputDir = $outputTo;
}
($fileName = $outputFile) =~ s/^.*?([^\/]+.pdf)$/\1/;

# Create data structure to read the results.
$dataSet->{'file'} = $galacticusFile;
$dataSet->{'store'} = 0;
&HDF5::Get_Parameters($dataSet);
&HDF5::Count_Trees($dataSet);
&HDF5::Select_Output($dataSet,0.0);

# Read galaxy data.
$dataSet->{'tree'} = "all";
&HDF5::Get_Dataset($dataSet,[
		       'volumeWeight',
		       'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega',
		       'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'
		   ]);
$dataSets  = $dataSet->{'dataSets'};

# Read the XML data file.
$xml     = new XML::Simple;
$data    = $xml->XMLin($galacticusPath."data/Morphological_Luminosity_Functions_2MASS_Devereux_2009.xml");

# Estimate bulge-to-total ratio ranges for each morphological class.
foreach $morphology ( @{$data->{'morphology'}} ) {
	# Get the luminosity function.
	$x  = pdl @{$morphology->{'magnitude'         }->{'datum'}};
	$y  = pdl @{$morphology->{'luminosityFunction'}->{'datum'}};
	$xMin    = where($x-5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'})-0.25                                                         ,$x == -23.25);
	$xMax    = where($x-5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'})+0.25                                                         ,$x == -23.25);
	$ySelect = where($y*($dataSet->{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'},$x == -23.25);
	$bulgeToTotal->{$morphology->{"class"}}->{"abundance"}              = $ySelect->index(0);
	$bulgeToTotal->{$morphology->{"class"}}->{"magnitude"}->{"minimum"} = $xMin   ->index(0);
	$bulgeToTotal->{$morphology->{"class"}}->{"magnitude"}->{"maximum"} = $xMax   ->index(0);
}

# Morphological classes (ordered).
@classes = ("Sc-Scd","Sb-Sbc","Sa-Sab","Lenticular","Elliptical");

# Structure giving bulge-to-total ranges for each morphological class.
$bulgeToTotal->{"Total"}->{"minimum"} = 0.00;
$bulgeToTotal->{"Total"}->{"maximum"} = 1.00;

# Compute cumulative morphological fractions.
$previousFraction = 0.0;
foreach $class ( @classes ) {
    unless ( $class eq "Total" ) {
	$bulgeToTotal->{$class}->{"abundance"} /= $bulgeToTotal->{"Total"}->{"abundance"};
	$bulgeToTotal->{$class}->{"abundance"} += $previousFraction;
	$previousFraction = $bulgeToTotal->{$class}->{"abundance"};
    }
}

# Compute cumulative fraction of model galaxies by bulge-to-total ratio.
$weight = where($dataSets->{'volumeWeight'},
		$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'} >= $bulgeToTotal->{"Total"}->{"magnitude"}->{"minimum"} &
		$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'} < $bulgeToTotal->{"Total"}->{"magnitude"}->{"maximum"}
		);
$ratio  = where($dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'},
		$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'} >= $bulgeToTotal->{"Total"}->{"magnitude"}->{"minimum"} &
		$dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'} < $bulgeToTotal->{"Total"}->{"magnitude"}->{"maximum"}
		);

$indices          = $ratio->qsorti;
$totalWeight      = $weight->sum;
$orderedRatios    = $ratio->index($indices);
$cumulativeWeight = ($weight->index($indices)->cumusumover)/$totalWeight;

# Interpolate to get corresponding model bulge-to-total ratios.
$previousFraction = 0.0;
foreach $class ( @classes ) {
    ($ratioBoundary,$error) = interpolate($bulgeToTotal->{$class}->{"abundance"},$cumulativeWeight,$orderedRatios);
    $bulgeToTotal->{$class}->{"minimum"} = $previousFraction;
    $bulgeToTotal->{$class}->{"maximum"} = $ratioBoundary;
    $previousFraction = $ratioBoundary;
}

# Open a pipe to GnuPlot.
my $plot;
my $gnuPlot;
my $plotFile = $outputFile;

# Loop over all morphologies and make plots.
my $iPlot = 0;
my @leafFiles;
my @plotFiles;
foreach $morphology ( @{$data->{'morphology'}} ) {

    if ( exists($bulgeToTotal->{$morphology->{'class'}}) ) {
	# Get the luminosity function.
	$x     = pdl @{$morphology->{'magnitude'              }->{'datum'}};
	$y     = pdl @{$morphology->{'luminosityFunction'     }->{'datum'}};
	$error = pdl @{$morphology->{'luminosityFunctionError'}->{'datum'}};
	
	# Scale for Hubble parameter.
	$x     += -5.0*log10($data->{'magnitudes'}->{'hubble'}/$dataSet->{'parameters'}->{'H_0'});
	$y     *= ($dataSet->{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'};
	$error *= ($dataSet->{'parameters'}->{'H_0'}/$morphology->{'luminosityFunction'}->{'hubble'})**$morphology->{'luminosityFunction'}->{'hubbleExponent'};

	# Construct Galacticus luminosity function.
	$magnitude = where($dataSets->{'magnitudeTotal:2MASS_Ks:observed:z0.0000:dustAtlas:vega'},
			     $dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'} >= $bulgeToTotal->{$morphology->{"class"}}->{"minimum"}
			   & $dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'} <= $bulgeToTotal->{$morphology->{"class"}}->{"maximum"}
	    );
	$weight    = where($dataSets->{'volumeWeight'},
			     $dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'} >= $bulgeToTotal->{$morphology->{"class"}}->{"minimum"}
			   & $dataSets->{'bulgeToTotalLuminosity:2MASS_Ks:observed:z0.0000:dustAtlas'} <= $bulgeToTotal->{$morphology->{"class"}}->{"maximum"}
	    );

	($yGalacticus,$errorGalacticus) = &Histograms::Histogram($x,$magnitude,$weight,differential => 1);
	
        # Compute chi^2.
	$thisChiSquared        = sum(($yGalacticus-$y)**2/($errorGalacticus**2+$error**2));
	$thisDegreesOfFreedom  = nelem($y);
	$chiSquared           += $thisChiSquared;
	$degreesOfFreedom     += $thisDegreesOfFreedom;

	# Make the plot.
	++$iPlot;
	(my $thisPlot = $plotFile) =~ s/\.pdf/$iPlot.pdf/;
	(my $thisPlotEPS = $thisPlot) =~ s/\.pdf$/.eps/;
	open($gnuPlot,"|gnuplot");
	print $gnuPlot "set terminal epslatex color colortext lw 2 solid 7\n";
	print $gnuPlot "set output '".$thisPlotEPS."'\n";
	print $gnuPlot "set xlabel '\$M_{\\rm K,vega}\$'\n";
	print $gnuPlot "set ylabel '\${\\rm d}n/{\\rm d}\log M_{\\rm K,vega} [\\hbox{Mpc}^{-3}]\$'\n";
	print $gnuPlot "set title '".$morphology->{'class'}." K Luminosity Function'\n";
	print $gnuPlot "set lmargin screen 0.15\n";
	print $gnuPlot "set rmargin screen 0.95\n";
	print $gnuPlot "set bmargin screen 0.15\n";
	print $gnuPlot "set tmargin screen 0.95\n";
	print $gnuPlot "set key spacing 1.2\n";
	print $gnuPlot "set key at screen 0.275,0.16\n";
	print $gnuPlot "set key left\n";
	print $gnuPlot "set key bottom\n";
	print $gnuPlot "set logscale y\n";
	print $gnuPlot "set mytics 10\n";
	print $gnuPlot "set format y '\$10^{\%L}\$'\n";
	print $gnuPlot "set xrange [-26.0:-19.0]\n";
	print $gnuPlot "set yrange [3.0e-6:3.0e-2]\n";
	print $gnuPlot "set pointsize 2.0\n";
	&PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $x,$y,
	    errorUp    => $error,
	    errorDown  => $error,
	    style      => "point",
	    symbol     => [6,7],
	    weight     => [5,3],
	    color      => $PrettyPlots::colorPairs{${$PrettyPlots::colorPairSequences{'slideSequence'}}[0]},
	    title      => $data->{'label'}.' [observed]'
	    );
	&PrettyPlots::Prepare_Dataset(
	    \$plot,
	    $x,$yGalacticus,
	    errorUp    => $errorGalacticus,
	    errorDown  => $errorGalacticus,
	    style      => "point",
	    symbol     => [6,7],
	    weight     => [5,3],
	    color      => $PrettyPlots::colorPairs{'redYellow'},
	    title      => 'Galacticus'
	    );
	&PrettyPlots::Plot_Datasets($gnuPlot,\$plot);

        # Close the pipe to GnuPlot.
	close($gnuPlot);
	&LaTeX::GnuPlot2PDF($thisPlotEPS);
	(my $leafName = $thisPlot) =~ s/^.*\/([^\/]+)$/$1/;
	push(@leafFiles,$leafName);
	push(@plotFiles,$thisPlot);
    }
}
&SystemRedirect::tofile("rm -f ".$outputFile."; cd ".$outputDir."; pdfmerge ".join(" ",@leafFiles)." tmp.pdf; cd -; mv ".$outputDir."/tmp.pdf ".$outputFile,"/dev/null");
&MetaData::Write($outputFile,$galacticusFile,$self);
unlink(@plotFiles);

# Output fit data.
if ( $showFit == 1 ) {
    $fitData{'name'} = "Cole et al. (2001) K-band luminosity function";
    $fitData{'chiSquared'} = $chiSquared;
    $fitData{'degreesOfFreedom'} = $degreesOfFreedom;
    $fitData{'fileName'} = $fileName;
    $xmlOutput = new XML::Simple (NoAttr=>1, RootName=>"galacticusFit");
    print $xmlOutput->XMLout(\%fitData);
}

exit;

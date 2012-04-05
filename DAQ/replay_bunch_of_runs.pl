: # feed this into perl *-*-perl-*-*
    eval 'exec perl -S $0 "$@"'
    if $running_under_some_shell;

use Getopt::Std;
use File::Basename;
use Env;

# Declare the subroutines
sub ltrim($);
sub rtrim($);


$executable = basename($0);

###  Get the option flags.
getopts('hvf:r:x:u:e:');

#GetOptions("h"=>\$opt_h,
#	   "v"=>\$opt_v,
#	   "f",\$opt_f,
#	   "r"=>\$opt_r);

###  Interpret the "help" flag
if ($opt_h==1){
    displayusage();
    exit;
}

###  Interpret the run number listing.
my @run_list = ();
$first_run = -1;
$last_run  = -1;
if ($opt_r ne ""){
    if ($opt_r =~ /,/){
	#  Comma seperated list.
	@run_list = map {int} sort {$a <=> $b} (grep {/^\d+$/ && $_>-1}
				      (split /,/, $opt_r));
	@discards = grep {!/^\d+$/ || $_<=-1 } (split /,/, $opt_r);
	if ($#discards > -1){
	    print STDERR
		"The following bad run numbers were discarded from the list: ",
		"@discards\n";
	}
    } elsif ($opt_r =~ /:/){
	# Colon seperated range.
	($first_run, $last_run) = split /:/, $opt_r, 2;
	if ($last_run <= $first_run){
	    print STDERR
		"The last run number is smaller than the first run number.",
		"  Discard the last run number\n";
	    $last_run = -1;
	    push @run_list, $first_run;
	} else {
	    #  Fill in the run list.
	    for ($i=$first_run; $i<=$last_run; $i++){
		push @run_list, $i;
	    }
	}
    } elsif ($opt_r =~ /^\d+$/) {
	#  Single run number.
	$first_run = $opt_r;
	$last_run = -1;
	push @run_list, $first_run;	
    } else {
	#  Unrecognisable value.
	print STDERR  "Cannot recognise the run number, $opt_r.\n";
	exit;
    }
} elsif($opt_f ne ""){
# -f flag, get the list of runs from command line
    $goodrunfile=$opt_f;
    if (-f $goodrunfile){
	print STDOUT "Using the run list file, $goodrunfile.\n";
    } else {
	print STDERR "Cannot find the good runs file \"$goodrunfile\".  Exiting.\n";
	exit;
    }
    @run_list = get_the_good_runs($goodrunfile);
} else {
    print STDERR  "ERROR: The \"run range\" must be specified.\n\n";
    displayusage();
    exit;
}

my $first_evt = -1;
my $last_evt = -1;
if ($opt_e =~ /:/){
    # Colon seperated range.                                                
    my ($tmp_first_evt, $tmp_last_evt) = split /:/, $opt_e, 2;
    if ($tmp_last_evt >= $tmp_first_evt && $tmp_first_evt>0){
	$first_evt = $tmp_first_evt;
	$last_evt = $tmp_last_evt;
    } 
}

foreach(@run_list){
    my $runno = $_;
    my $runlimit5 = 10000;
    my $command;
    my $hbookout;
    my $rootout;
    my $rawfile;
    my $zipfile;
    my $odbfile;
    my $runno5;
    

    printf("\"$runno\"\n");
    
    if ($runno < $runlimit5){
	$runno5 = "0$runno";
    } else { #10000 or more
	$runno5 = $runno;
    }
    
    $inputdatadir = "/home/data_raw/2011";
    $outputhbookdir = "/home/data_analyzed/2011/hbooks";
    $outputrootdir = "/home/data_analyzed/2011/rootfiles";
    $hbookout = "$outputhbookdir/full$runno5.rz";
    $rootout = "$outputrootdir/full$runno5.root";
    $rawfile = "$inputdatadir/run$runno5.mid";
    $zipfile = "$inputdatadir/run$runno5.mid.gz";
    
    if($opt_u ne "") {
	$odbfile = $opt_u;
    } else {
	$odbfile = "$inputdatadir/run$runno5.odb";
    }
    
    #now loop for raw data file
    if (-e $rawfile)
    { 
	#good, do nothing
    } elsif (-e $zipfile) {
	printf("Raw file is zipped into $zipfile!!!!\n");
	printf("Unzipping it .... Please be patient\n");
	system("gunzip $zipfile");
    } else {
	printf("Neither $rawfile nor $zipfile exist!!!! Exit now!!!!\n");
	next;
    }
    
    # split on different version of analyzer
    my $executable;
    if($opt_x eq "on") {
	$executable = "analyzer_Before2008Sep23On";
    } elsif ($opt_x eq "off") {
	$executable = "analyzer_Before2008Sep23Off";
    } else {
	$executable = "analyzer";
    }

    if($first_evt!=-1 && $last_evt!=-1) {
	$command = "./$executable -n $first_evt $last_evt -w  -i $rawfile -o $hbookout; h2root $hbookout $rootout";
    } else {
	$command = "./$executable -w  -i $rawfile -o $hbookout; h2root $hbookout $rootout";
    }
    print("$command\n");
    if($opt_v==1) {
	print("$command\n");
    }

    # search for if there is another analysis process which is running
    # if so, do nothing
    my $grepcmd = "ps -ef | grep -i analyzer | grep $runno | grep -v grep";
    chomp(my $output = qx!$grepcmd!);
    print "####### $grepcmd\n";
    print "####### $output\n"; 

    my $offline_dir = "/home/daq/offline";
    
    if ($output) {
	print "An analyzer has already been running on $runno. Do nothing!\n";
    }  else {
	#create a temporary dir for .SHM file
	my $tmpdir = "$offline_dir/tmp$runno";
	system("mkdir $tmpdir");
	$ENV{MIDAS_DIR} = $tmpdir;
	
	#modify LD_LIBRARY_PATH
	my $oldenv = $ENV{LD_LIBRARY_PATH};
	$ENV{LD_LIBRARY_PATH}="$oldenv:$offline_dir:";
	my $newenv = $ENV{LD_LIBRARY_PATH};
	print "$newenv\n";
	
	#load the correct ODB file
	system("odbedit -c \"load $odbfile\"");
	print "$odbfile Loaded!!!!!\n";
	
	#replay the run
	system("$command");
	printf("gzipping $hbookout ... \n");
	system("gzip -f $hbookout");
	system("gzip -f $rawfile");
	
	system("rm -rf $tmpdir");
    }
}

exit;

################################################
sub displayusage {
    print STDERR
	"\n",
	"$executable provides a means of multi-run OFFLINE\n",
	"replay given a list of run numbers. It also does the h2root conversion.\n",
	"Usage:\n\t$executable [-h]\n",
	"\t$executable [-v] [-f <file of run list>] [-u <odbfile>] [-x <on/off>]\n",
	"\t$executable [-v]  [-r <run range>] [-u <odbfile>] [-x <on/off>]\n",
	"Options:\n\t-h\n\t\tPrint usage information\n",
	"\t-v\n",
	"\t\tEnable verbose output.\n",
	"\t-f <file of run list>\n",
	"\t\tThis specifies the name of  the file containing the\n",
	"\t\tlist of runs to be replayed. If the full path is not given\n",
	"\t\tthe file will be searched for in the same directory\n",
	"\t\tthat contains  the executable.It is a plain text file\n",
	"\t\tcontaining one run number per line.\n",
	"\t\tSee an example file at \"runs2replay.lst\".\n",
	"\t-r <run range>\n",
	"\t\tThis flag  specifies the numbers  of the runs to be\n",
	"\t\tstaged in from tape.  Three formats are permitted:\n",
	"\t\t  A single run number:      15789\n",
	"\t\t  A comma separated list:   15775,15781,15789\n",
	"\t\t  A colon separated range:  15775:15793\n",
	"\t-x <on/off>\n",
	"\t\t Use use supplied flipper state to blind the clock\n",
	"\t-u <odbfile>\n",
	"\t\tUse user supplied ODB file\n",
	"\t-e <first_evt:last_evt>\n",
        "\t\tUse colon to separate event range:  1:10000\n";
}

sub get_the_good_runs ($) {
    local(@filelist, $goodfile, $runnumber);
    ($goodfile) = @_;
    @filelist = ();

    open GOODRUNS, $goodfile;
    %good_run = ();
    while (<GOODRUNS>){
	next if /^#/;
	next if /^$/;
	chomp;
	$runnumber = $_;
	$runnumber = ltrim($runnumber);
	$runnumber = rtrim($runnumber);
	push @filelist, $runnumber;
    }
    return @filelist;
}

# Left trim function to remove leading whitespace
sub ltrim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    return $string;
}
# Right trim function to remove trailing whitespace
sub rtrim($)
{
    my $string = shift;
    $string =~ s/\s+$//;
    return $string;
}

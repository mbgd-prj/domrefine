#!/usr/bin/perl -w
use strict;
use File::Basename;
use Getopt::Std;
my $PROGRAM = basename $0;
my $USAGE=
"Usage: $PROGRAM 'COMMAND LINE'
OPTION:
 -q QUEUE
 -N JOB_NAME
 -l RESOURCE
DOMREFINE_QUEUE is used as the default queue when -q is not specified
DOMREFINE_QSUB_TMP is used as SGE in/out directory
DOMREFINE_QSUB_OPT is used as additional qsub options
";
# DOMREFINE_PROJECT_NAME is used for qsub -P option

use DomRefine::Read;

my %OPT;
getopts('q:N:l:', \%OPT);

@ARGV == 0 and die $USAGE;

my $JOB_NAME = $OPT{N} || $PROGRAM;

### Create temporary script ###
my $TMP_SCRIPT = "$ENV{DOMREFINE_QSUB_TMP}/$ENV{HOSTNAME}.$$.$PROGRAM.$JOB_NAME.sh";
open(SCRIPT, ">$TMP_SCRIPT") || die "Cannot open $TMP_SCRIPT";
print SCRIPT '#$ -S /bin/sh', "\n";
print SCRIPT "@ARGV\n";
close(SCRIPT);

### Submit job ###
my $QSUB = "qsub";
$QSUB .= " -e $ENV{DOMREFINE_QSUB_TMP} -o $ENV{DOMREFINE_QSUB_TMP}";
$QSUB .= " -V -N $JOB_NAME $ENV{DOMREFINE_QSUB_OPT}";

if ($OPT{q}) {
    $QSUB .= " -q $OPT{q}";
} elsif ($ENV{DOMREFINE_QUEUE}) {
    $QSUB .= " -q $ENV{DOMREFINE_QUEUE}";
}

if ($OPT{l}) {
    $QSUB .= " -l $OPT{l}";
}

system "$QSUB $TMP_SCRIPT";

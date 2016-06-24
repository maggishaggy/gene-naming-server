#!/usr/local/bin/perl -w

print "Content-type:  text/html\n\n";

use CGI qw(:standard);

print "Testing Remote Submission<BR>";

print "<BR>";
print "<BR>";

my $BSUB_COMMAND='hostname';

my $COMMAND = ". /opt/lsf/conf/profile.lsf\; bsub -u bob_freeman\@hms.harvard.edu -W 5 \'" . $BSUB_COMMAND . "\'  ";


my $ret = system('whoami');
print "COMMAND: <BR>";
print $COMMAND . "<BR>";


print "<BR>";
print "<BR>";

my $RESULT = qx/$COMMAND/;

print "COMMAND RESULT = " .$RESULT . "<BR>" ;


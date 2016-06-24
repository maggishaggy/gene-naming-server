#!/usr/bin/perl
#
#
# CGI script called from make_mz_ref_db web page. Will execute:
#   1. uploads file, verifying size and # transcripts
#   2. file is copied (if not already) to scratch for work
#   3. make protein ref scripts are called...
# 
  

# DETAILS
#   get new, unique directory on scratch
#   grab upload filename and transfer to unique name on scratch
#   verify filesize and throw error if too large
#   verify # transcripts and throw error if too large
#   start mz_ref script part1
#
# 11/13/13, rmf
#   added email address validator (well-formed)
#
# 11/15/13, rmf
#   sends out email after LSF job submitted
#   generates stats of input file
#   hangs on large files; move inline analysis to LSF mini immediate jobs
#
# 11/19/13, rmf
#   better error handling vis-a-vis browser
#
# 11/21/13, rmf
#   fixed job log not logging emails properly
#   fixed success emails not having proper format
#   moved counting sequences out of cgi script and moved to www_make_mz_protein_ref_part1
#   changed html message to reflect updated email message to be sent
#
# 12/2/13, rmf
#   added routine to clean up temp space (our 'local' 'scratch' space)
#
#   12/16/13, rmf
#       added options for clust_pcnt_ID, clust_overlap_len, and asm_pcnt_ID
#       included options in outgoing email
#       added extra field to command log to track if submission coming from DEV or PROD

use strict;

use lib "$ENV{'DOCUMENT_ROOT'}/../support/mz_ref_db/lib/perl5";

use English '-no_match_vars';
use FindBin;
use File::Temp qw/ tempfile tempdir /;
use CGI qw(:standard escape escapeHTML);
use File::Basename;
use Data::Dumper;
#use lib "/www/kirschner.med.harvard.edu/support/home/local";
#use lib '/groups/kirschner/local/lib';
use Email::Valid;
use Email::Simple;
use Email::Sender::Simple qw(sendmail);

my $PROJECT_HOME = "$ENV{'DOCUMENT_ROOT'}/../support/mz_ref_db";
$ENV{'PERL5LIB'} = "$PROJECT_HOME/lib/perl5";
$ENV{'PROJECT_HOME'} = $PROJECT_HOME;


my $debug = 0;
my ($cmd, $result);
my $user_fasta_filename;
my ($clust_pcnt_ID, $clust_min_ovrlap_len, $asm_pcnt_ID);
my $email_body;
#my $BIN_ROOT = "$FindBin::Bin/..";
my $BIN_ROOT = "$PROJECT_HOME/bin";
#my $KLAB_BIN='/groups/kirschner/local/bin';
my $max_file_size = 500 * 1024 * 1024; # 500 MB limit
$CGI::POST_MAX = $max_file_size;

my $dir_template = 'tempdirXXXXXXXX';
my $file_template = 'tempfileXXXXXXXX';
my $scratch = '/groups/kirschner_www/mz_ref_db'; 
my $job_log_file = 'job_log.txt';
my $blast_db;           # place in support/mz_ref_db/lib folder of this module
my $tln_E_filter;       # E cutoff for translation filter

my ($running_lsf_jobs, $total_lsf_jobs);
 
# All output going to user's web browser
my $browser = CGI->new;

print header('text/html'),
      start_html('Making a Mass Spec Protein Reference...');
#print "Content-type: text/html\n\n";
#print "<HTML><BODY>";
  
if (0) {
    print "%ENV:<BR><BLOCKQUOTE>";
    print Dumper(%ENV);
    print "</BLOCKQUOTE><BR><BR>";
}

# grab form params, or supply our own if debugging
my ($UPLOAD_FH, $sender_email);
if (!$debug) {
    #$UPLOAD_FH             = param('upload_file');
    $UPLOAD_FH              = upload('upload_file');
    $user_fasta_filename    = param('upload_file');
    $sender_email           = param('email');
    $clust_pcnt_ID          = param('ClustPcntID');
    $clust_min_ovrlap_len   = param('ClustMinOvrlapLen');
    $asm_pcnt_ID            = param('AsmPcntID');
    $blast_db               = param('blast_db');
    $tln_E_filter           = param('TlnEFilter');
} else {
    $UPLOAD_FH = 'XGI.shuff.small.fa';
    $sender_email = 'bob_freeman@hms.harvard.edu';
    $clust_pcnt_ID = 94;
    $clust_min_ovrlap_len = 100;
    $asm_pcnt_ID = 93;
}

#$result = `date`;
#print "$result ...<BR><BR>";
print `date` . "<BR><BR>";

# print "(Please be patient... DO NOT HIT REFRESH button or else your job submission will be lost!)<BR><BR>$result ...<BR><BR>";
# print "Validating email address...";
# first check to see that email is well-formed
if (! Email::Valid->address($sender_email)) {
    # The email address is not valid
    exit_with_error_msg('I\'m sorry, but your e-mail address does not appear to be valid. Please use your browser\'s Back button and correct the entry.');
}


#print "Validating FASTA file data...";
# have a supplied filename. let's do our magic!
# create the temp dir and temp file

if (!defined $UPLOAD_FH) {
    exit_with_error_msg('I\m sorry, but the FASTA filename seems to be missing. Please use your browser\'s Back button and specify a valid filename.');
}


# make our scratch dir if not already there
if (! -e $scratch) {
    if (! mkdir $scratch) {
        exit_with_error_msg("I'm sorry, but there was a fatal error creating the scratch directory $scratch:<BR>"
            . "msg: $!<BR> Please contact support with this error message.");
    }
}

# create tempdir and tempfile, opening the latter when created
my ($tempdir_path, $tempdir, $DATFH, $fasta_pathname, $working_dir_path, $fasta_filename);
eval {
    $tempdir_path = tempdir( $dir_template, DIR => $scratch, CLEANUP => 0 );
    `chmod go+rx $tempdir_path`;
};
exit_with_error_msg("I'm sorry, but there was a fatal error creating the temporary for processing. "
    . "Please contact support with this error message.") if $EVAL_ERROR;
$tempdir = basename($tempdir_path);
$working_dir_path = $tempdir_path;

eval {
    ($DATFH, $fasta_pathname) = tempfile( $file_template, SUFFIX => '.fasta', DIR => $tempdir_path, UNLINK => 0 );
};
exit_with_error_msg("I'm sorry, but there was a fatal error creating the scratch directory $scratch:<BR>"
    . "Please contact support with this error message.") if $EVAL_ERROR;
($fasta_filename) = fileparse($fasta_pathname);


# and now  copy the file up
my ($data, $length, $chunk);
if ($debug) {
    my $DEBUGFH;
    if (! open ($DEBUGFH, "<$UPLOAD_FH")) {
        exit_with_error_msg("Error opening file $UPLOAD_FH for reading: $!");
    }
    $UPLOAD_FH = \*$DEBUGFH;
}
while ($chunk = read ($UPLOAD_FH, $data, 1024)) {
    print $DATFH $data;
    
    $length += $chunk;
    if ($length > $max_file_size) {
        close $DATFH;
        unlink $fasta_pathname;
        exit_with_error_msg('I\'m sorry, but the upload file is too large. The file size limit is 500 MB. '
                             .' Please use your browser\'s Button and select a smaller file.');
    }
}
close $DATFH;

if (! chdir $working_dir_path) {
    exit_with_error_msg("I'm sorry but there was an error changing to the working dircectory $working_dir_path:<BR>"
        . "msg: $!<BR> Please contact support with this error message.");
}


# now sanity check # of sequences.
# park this in www_make_mz_protein_ref_part1 script, as it's taking too long to count them...
my $num_seqs = 0;
# $num_seqs = `grep -c ">" $fasta_filename`;
# $num_seqs += 0; 
# if ($num_seqs > $max_num_seqs) {
#     unlink $fasta_pathname;
#     exit_with_error_msg('I\'m sorry, but your file contains more than the allowed number of sequences (500K). '
#                              .' Please adjust your file, use your browser\'s Button and submit again.');
# }


# and call make_mz_ref_script....
#   non-server version of this script is make_mz_protein_ref_from_mRNA.pl
$cmd = <<"BSUB_CMD";
 . /opt/lsf/conf/profile.lsf\; bsub -q sysbio_7d -J www_mz_ref_db \\
     -n 8 -R "span[hosts=1]" \\
     -N -o $fasta_filename.stdout \\
     -e $fasta_filename.stderr \\
 $BIN_ROOT/mz_ref_db/www_make_mz_protein_ref_part1.pl -c 8 -i $fasta_filename -o final_proteins.fa \\
     -p $clust_pcnt_ID -v $clust_min_ovrlap_len -a $asm_pcnt_ID \\
     -b $blast_db --email $sender_email -f $tln_E_filter
BSUB_CMD
;
chomp $cmd;                             #(remove the extra line ending so that we tack on extra params...
my ($lsf_sub_result, $lsf_job_ID);
$lsf_sub_result	= `$cmd`;               # Job <1636426> is submitted to queue <mini>.
($lsf_job_ID) = $lsf_sub_result =~ m/Job <(\d+)>/;
if (defined $lsf_job_ID) {

    # save job ID in work directory
    `touch jobID.$lsf_job_ID.jobID`;
    # modify job name to include Job ID on the end
    #`. /opt/lsf/conf/profile.lsf\; bmodify -J www_mz_ref_db_$lsf_job_ID $lsf_job_ID`;

    # get LSF status
    ($running_lsf_jobs, $total_lsf_jobs) = current_lsf_job_status(); 
    
    my ($success_html, $success_email) = init_heredoc();
    print $success_html;

    # now send out same message as email
    my $email = Email::Simple->create(
        header => [
          To      => "$sender_email",
          From    => 'Bob Freeman <bob_freeman@hms.harvard.edu>',
          Subject => "Mass spec protein reference creation jobID $lsf_job_ID started",
        ],
        body => $success_email,
    );
    # now send the message
    if (0) {
        sendmail($email);
    }

    # debug
    # save email as text so we know what the problem is...
    open (my $sender_emailFH, ">email.txt");
    print {$sender_emailFH} $success_email;
    close $sender_emailFH;
    # debug
    #   and do the same thing for the job command!
    open (my $commandFH, ">command.txt");
    print {$commandFH} $cmd;
    close $commandFH;
    

    # log job entry
    $result = `date`;
    chomp $result;
    my $website = 'prod';
    $website = 'dev' if $ENV{'DOCUMENT_ROOT'} =~ m|/dev|; 
    open (my $LOGFH, ">>$scratch/$job_log_file");
    print {$LOGFH} join("\t", $result,$lsf_job_ID,$tempdir,$fasta_filename,$num_seqs,$sender_email,$ENV{'REMOTE_ADDR'},$website) . "\n";
    close $LOGFH;


} else {
    exit_with_error_msg("I'm sorry, but there was a fatal error submitting the job to LSF.<BR>"
            . "msg: $lsf_sub_result<BR>Please contact support with this error message.");
}

# 12/2/13, rmf
# added cleanup routine for working directory
`. /opt/lsf/conf/profile.lsf\; bsub -q mini "find $scratch -mtime +7 -delete"`;
#


exit; 

#
#
# END OF PROGRAM

sub exit_with_error_msg {

    my ($message) = @_;
    print $message, end_html;
    exit;
}

sub init_heredoc {

    my $success_email = <<"_SUCCESS_";
Success! Your job has been submited to the pipeline for processing with the following options:

Filename:           $user_fasta_filename
# sequences:        TBD

Cluster % ID:       $clust_pcnt_ID
Clust overlap len   $clust_min_ovrlap_len
Assembly % ID:      $asm_pcnt_ID

BLAST database:     $blast_db
Translate E cutoff:  $tln_E_filter

JobID:              $lsf_job_ID
Contact e-mail:     $sender_email
Job log:            http://kirschner.med.harvard.edu/downloads/mz_ref_db_results/$tempdir/$fasta_filename.stdout
Error log:          http://kirschner.med.harvard.edu/downloads/mz_ref_db_results/$tempdir/$fasta_filename.stderr

Cmd:                $cmd


An updated copy of this message will be sent via email and will include both validation results and an estimate of the compute time required. Please be aware that currently there are $total_lsf_jobs jobs ahead of you, $running_lsf_jobs of which are running.  Any questions, please do not hesitate to contact us at bob_freeman\@hms.harvard.edu.
_SUCCESS_
;

    my $success_html = $success_email;
    $success_html =~ s/\n/<BR>\n/g;
    return ($success_html, $success_email);

} 

sub current_lsf_job_status {

    # calculate how many jobs are in the queue total,
    my $total_lsf_jobs = `. /opt/lsf/conf/profile.lsf\; bjobs -w -q sysbio_7d -J www_mz_ref_db | grep www_mz_ref_db | wc -l`;
    ($total_lsf_jobs) = $total_lsf_jobs =~ m/^(\d+)/;
    #print STDERR "***** LSF total $total_lsf_jobs\n";

    # and now how many running
    my $running_lsf_jobs = `. /opt/lsf/conf/profile.lsf\; bjobs -w -q sysbio_7d -J www_mz_ref_db | grep RUN | wc -l`;
    ($running_lsf_jobs) = $running_lsf_jobs =~ m/^(\d+)/;
    #print STDERR "***** LSF running $running_lsf_jobs\n";

    return ($running_lsf_jobs, $total_lsf_jobs);
}
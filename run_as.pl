use strict;
use warnings;
use Cwd;

use File::Basename;
use lib dirname (__FILE__);

use Samples;
use Analysis;

## To run AS analysis
## Files needed:
## - Sample Mapping
## - Sample Grouping
## - Sample Comparison
## - Alignments folder (where STAR output is)
## - Head output dir
## - Genome location (fastq.gz) (defaulted)
## - GTF location (gtf and indexed gtf) (defaulted)
## - SpliceTools annotation location (defaulted)

use Getopt::Long qw( GetOptions );

# Set defaults in args:
my %args = ( 
    "ref_fasta" => "abc",
    "gtf_file" => "def",
    "gtf_index_file" => "zyx",
    "anno_file" => "ghi",
    "output_dir" => getcwd
);

%args = get_args(%args);
check_required_args(%args);
my %samp = Samples::setup_sample_hash(
    $args{mapping_file},
    $args{grouping_file},
    $args{comparison_file},
    $args{alignments_dir}
);

my $read_len = get_read_len(\%samp);

Analysis::run_kallisto(\%args, \%samp);
Analysis::run_rmats(\%args, \%samp, $read_len);
Analysis::run_splicetools(\%args, \%samp);


sub get_args {
    my %args = @_;
    my $help = 0;
    GetOptions( \%args,
        'mapping_file=s',
        'grouping_file=s',
        'comparison_file=s',
        'alignments_dir=s',
        'ref_fasta=s',
        'gtf_file=s',
        'gtf_index_file=s',
        'annotation_file=s',
        'strand=s',
        'anno_file=s',
        'output_dir=s'
    ) or die "something went wrong";

    return %args;
}

sub check_required_args {
    my %args = @_;
    my @missing_req_args;

    foreach my $req_arg ("mapping_file", "grouping_file", "comparison_file", "alignments_dir", "strand", "annotation_file"){
        if(! exists $args{$req_arg}){
            push(@missing_req_args, $req_arg)
        } else {
            if($req_arg =~ /file/){
                if(! -e $args{$req_arg}){
                    die $req_arg . " is not a file: " . $args{$req_arg};
                }
            } elsif($req_arg =~ /dir/) {
                if(! -d $args{$req_arg}){
                    die $req_arg . " is not a directory: " . $args{$req_arg};
                }
            }
        }
    }

    if(@missing_req_args){
            die "Missing required arguments: " . join(", ", @missing_req_args);
        }
}

# Assuming reads are in fastq.gz format.
# Also assuming every sample has the same read_len
sub get_read_len {
    my $samp = shift;
    
    my $fq = glob $samp->{ $samp->{$samp_list}->[0] }->{fq} . "/*fastq.gz";

    my $fq_line = `zcat < $fq | sed '2q;d'`;
    chomp $fq_line;
    
    return length($fq_line);
}



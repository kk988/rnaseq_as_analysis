package Samples;

use strict;
use warnings;

sub setup_sample_hash {
    my $mapping = shift;
    my $grouping = shift;
    my $comparisons = shift;
    my $align_dir = shift;

    my %samp;

    # Grab sample name, and mapping dirs - set into samp{samp_nam}{fq}
    process_mapping_file(\%samp, $mapping);
    process_grouping_file(\%samp, $grouping);
    process_comparison_file(\%samp, $comparisons);
    assign_alignments(\%samp, $align_dir);

    return %samp
}

sub process_mapping_file {
    my $samp = shift;
    $samp->{samp_list} = [];
    my $mapping = shift;

    open(MF, '<', $mapping) or die $!;
    foreach my $line (<MF>){
        chomp $line;
        my (undef, $sname, undef, $fq_dir,$ended) = split("\t", $line);

        if($sname ~~ @{ $samp->{samp_list} }){
            if( $samp->{$sname}{ended} ne $ended ){
                die "Sample $sname reads are both SE and PE - remove one to create consistency";
            }

            push(@{ $samp->{$sname}->{fq} }, $fq_dir);
        } else { 
            push( @{$samp->{samp_list}}, $sname);
            $samp->{$sname} = {
                'fq' => [$fq_dir],
                'ended' => $ended
                };
        }
    }
    close(MF);
}

sub process_grouping_file {
    my $samp = shift;
    my $grouping = shift;

    $samp->{group} = {};

    open(GF, '<', $grouping) or die $!;
    foreach my $line (<GF>){
        chomp $line;
        my ($sname, $group) = split("\t", $line);

        if( ! exists $samp->{group}->{$group}){
            $samp->{group}->{$group} = [$sname];
        } else {
            push( @{$samp->{group}->{$group}}, $sname );
        }
    }
    close(GF);
}

sub process_comparison_file {
    my $samp = shift;
    my $comparison = shift;

    $samp->{comps} = {};

    open(CF, '<', $comparison) or die $!;
    foreach my $line (<CF>){
        chomp $line;
        $line =~ s/\r//g;
        
        my ($c1, $c2) = split("\t", $line);
        
        if( ! exists $samp->{comps}->{$c1}){
            $samp->{comps}->{$c1} = [$c2];
        } else {
            push( @{ $samp->{comps}->{$c1} }, $c2);
        }
    }
    close(CF);
}

sub assign_alignments {
    my $samp = shift;
    my $align_dir = shift;

    foreach my $sname (@{ $samp->{samp_list} }){
        my @bams = glob("$align_dir/*$sname.bam");

        if( ! @bams){
            print "Sample $sname is missing bams!";
        }

        $samp->{$sname}->{bams} = @bams; 
    }
}

1;
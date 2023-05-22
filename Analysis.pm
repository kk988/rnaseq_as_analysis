package Analysis;

use strict;
use warnings;
use File::Basename;

sub run_kallisto {
    my $args = shift;
    my $samp = shift;

    # for each sample, run kallisto
    mkdir $args->{output_dir} . "/kallisto";
    foreach my $sample (@{$samp->{samp_list}}){
        print "Running Kallisto for sample $sample\n";

        my $k_out = $args->{output_dir} . "/kallisto/" . $sample;
        $samp->{$sample}->{kallisto_dir} = $k_out;
        mkdir $k_out;

        my $fqs = grab_fq_pairs($samp->{$sample}->{fq}, $samp->{$sample}->{ended});

        my $stranded = "";
        if($args->{strand} =~ /reverse/i ){
            $stranded = " --rf-stranded";
        } elsif ($args->{strand} =~ /forward/i ){
            $stranded = " --fr-stranded";
        }

        my $cmd = "singularity exec --bind " . join(" --bind ", @{$samp->{$sample}->{fq}}) 
            . " --bind " . $args->{output_dir}
            . " docker://zlskidmore/kallisto:0.48.0 /usr/local/bin/kallisto quant"
            . " -i " . $args->{gtf_index_file}
            . " -o " . $k_out
            . " -b 100 -t 2" . $stranded . $fqs;

        print $cmd . "\n";

    }
}

# Some assumptions here:
# we will always deal with fastq.gz
# They have _R1_, and the samples don't have that pattern
sub grab_fq_pairs {
    my $fq_dir = shift;
    my $pairedness = shift;

    my @r1_reads = glob("$fq_dir/*_R1_*fastq.gz");

    if( $pairedness ne 'PE'){
        return join(" ", @r1_reads);
    }

    my $fq_str = "";
    foreach my $read (@r1_reads) {
        my $r2_read = $read =~ s/_R1_/_R2_/;
        $fq_str .= " $read $r2_read";
    }

    return $fq_str;
}

# rmats runs by comparison. 
# You can put "paired" analysis meaning that
# The samples in b1 are paired with the cooresponding
# b2 samples. Because we don't have a way to pair
# the samples, we aren't using paired option.
sub run_rmats {
    my $args = shift;
    my $samp = shift;
    my $read_len = shift;

    my $mats_out = $args->{output_dir} . "/mats";
    my $mats_tmp = "$mats_out/tmp";

    my $stranded = "";
        if($args->{strand} =~ /reverse/i ){
            $stranded = " --libType fr-secondstrand";
        } elsif ($args->{strand} =~ /forward/i ){
            $stranded = " --libType fr-firststrand";
        }

    # Should we add --variable-read-length ??
    my $cmd =  "singularity exec --bind " . $args->{alignments_dir} 
        . " --bind " . dirname $args->{gtf_file} 
        . " --bind $mats_out docker://xinglab/rmats:v4.1.2 python /rmats/rmats.py"
        . " --gtf " . $args->{gtf_file} . " $stranded"
        . " --readLength $read_len "
        . " --nthread 4"
        . " --tmp $mats_tmp";

    foreach my $c1 keys %{$samp->{comps}} {
        foreach my $c2 @{$samp->{comps}->{$c1}} {
            my $result_dir = "$mats_out/$c1\_vs_$c2";
            mkdir $result_dir;
            my $b1_file = "$result_dir/b1.txt";
            my $b2_file = "$result_dir/b2.txt";

            my $b1_content = join(",", find_bams($args->{alignments_dir}, $samp->{group}->{$c1}));
            open(BF, '>', $b1_file);
            print BF $b1_content;
            close(BF);

            my $b2_content = join(",", find_bams($args->{alignments_dir}, $samp->{group}->{$c2}));
            open(BF, '>', $b2_file);
            print BF $b2_content;
            close(BF);

            $s_cmd = $cmd . " --b1 $b1_file --b2 $b2_file --od $result_dir";

            print "$s_cmd\n";
        }
    }

}

sub find_bams {
    my $align_dir = shift;
    my $sample_ary = shift;

    my @bams;
    
    foreach my $samp @{$sample_ary}{
        my @sample_bams = glob "$align_dir/*_$samp.bam";
        if(scalar @sample_bams != 1){
            die "Unable to find correct number of bams for $samp: " . join(",", @sample_bams);
        }
        push(@bams, $sample_bams[0]);
    }

    return @bams;
}

sub run_splicetools {
    my $args = shift;
    my $samp = shift;

    my $cmd = "singularity exec --bind " . $args->{output_dir}
        . " --bind " . dirname $args->{ref_fasta} 
        . " --bind " . dirname $args->{annotation_file}
        . " docker://kk988/splicetools:v1.1_1 perl /usr/bin/SpliceTools-1.1/bin/RIMedley.pl"

    foreach my $c1 keys %{$samp->{comps}} {
        foreach my $c2 @{$samp->{comps}->{$c1}} {
            my $mats_dir = "$mats_out/$c1\_vs_$c2";
            my $expr_file = combine_expression_files($samp, $c1, $c2, $mats_dir);
            my $ri_file = "$mats_dir/RI.MATS.JCEC.txt";
            my $se_file = "$mats_dir/SE.MATS.JCEC.txt";
            my $num_samples = join(",", scalar(@{$samp->{group}->{$c1}}), scalar(@{$samp->{group}->{$c2}}));

            my $ri_medley = $cmd 
                . " -r $ri_file"
                . " -a " . $args->{annotation_file}
                . " -g " . $args->{ref_fasta}
                . " -e $expr_file"
                . " -TPM 2,2"
                . " -SN $num_samples"
                . " -f 0.01";

            print "$ri_medley\n";

            my $se_medley = $cmd 
                . " -s $se_file"
                . " -a " . $args->{annotation_file}
                . " -g " . $args->{ref_fasta}
                . " -e $expr_file"
                . " -TPM 2,2"
                . " -SN $num_samples"
                . " -f 0.01";

            print $se_medley;

        }
    }
}


# It's easier to do bash manipulation of files than
# use perl, and it should use less memory.
sub combine_expression_files {
    my $samp = shift;
    my $grp1 = shift;
    my $grp2 = shift;
    my $out_dir = shift;

    my $expr = "$out_dir/expr.tsv";
    my $tmp_expr = "$out_dir/tmp_expr.tsv";
    my $final_expr = "$out_dir/expression_file.tsv"
    my $header_line = "target_id";

    my $first_samp = $samp->{samp_list}->[0];
    # First grab column
    my $cmd = "cut -f 1 " . $samp->{ $first_samp }->{kallisto_dir} . "/abundance.tsv > $expr";
    `$cmd`;

    # Then iteratively grab the column from each sample's abundance file
    foreach my $sname (@{$samp->{group}->{$grp1}}){
        my $abundance = $samp->{ $sname }->{kallisto_dir} . "/abundance.tsv";
        $header_line .= " $sname";
        my $cmd2 = "awk '{print \$NF}' $abundance | paste $expr - > $tmp_expr; mv $tmp_expr $expr;";
        `$cmd`;
    }

    foreach my $sname (@{$samp->{group}->{$grp2}}){
        my $abundance = $samp->{ $sname }->{kallisto_dir} . "/abundance.tsv";
        $header_line .= " $sname";
        my $cmd2 = "awk '{print \$NF}' $abundance | paste $expr - > $tmp_expr; mv $tmp_expr $expr;";
        `$cmd`;
    }

    # Finally - add the header line:
    my $cmd3 = "echo \"$header_line\" | tr ' ' '\t' > $final_expr; tail -n +2 $expr >> $final_expr;";

    return $final_expr;
}

1;
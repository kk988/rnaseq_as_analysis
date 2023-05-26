package Analysis;

use strict;
use warnings;
use File::Basename;

sub run_kallisto {
    my $args = shift;
    my $samp = shift;

    # for each sample, run kallisto
    mkdir $args->{output_dir} . "/kallisto";

    my $log_dir = $args->{output_dir} . "/logs";
    mkdir $log_dir;

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

        my $cmd = "singularity exec --bind " . join(" --bind ", @{$samp->{$sample}->{fq}});
        $cmd .= " --bind " . $args->{output_dir};
        $cmd .= " --bind " . dirname $args->{gtf_index_file};
        $cmd .= " docker://zlskidmore/kallisto:0.48.0 /usr/local/bin/kallisto quant";
        $cmd .= " -i " . $args->{gtf_index_file};
        $cmd .= " -o " . $k_out;
        $cmd .= " -b 100 -t 2 $stranded $fqs";
        $cmd .= " > $log_dir/kallisto_$sample.log 2>&1";

        print $cmd . "\n\n";
        system($cmd);

    }
}

# Some assumptions here:
# we will always deal with fastq.gz
# They have _R1_, and the samples don't have that pattern
sub grab_fq_pairs {
    my $fq_dirs = shift;
    my $pairedness = shift;

    my $fq_str = "";

    foreach my $fq_dir (@{$fq_dirs}){
        my @r1_reads = glob("$fq_dir/*_R1_*fastq.gz");

        if( $pairedness ne 'PE'){
            return join(" ", @r1_reads);
        }


        foreach my $read (@r1_reads) {
            (my $r2_read = $read) =~ s/_R1_/_R2_/;
            $fq_str .= " $read $r2_read";
        }
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
    mkdir $mats_out;
    my $mats_tmp = "$mats_out/tmp";
    mkdir $mats_tmp;

    my $log_dir = $args->{output_dir} . "/logs";
    mkdir $log_dir;

    my $stranded = "";
        if($args->{strand} =~ /reverse/i ){
            $stranded = " --libType fr-secondstrand";
        } elsif ($args->{strand} =~ /forward/i ){
            $stranded = " --libType fr-firststrand";
        }

    # Should we add --variable-read-length ??
    my $cmd =  "singularity exec --bind " . $args->{alignments_dir}; 
    $cmd .= " --bind " . dirname $args->{gtf_file} ;
    $cmd .= " --bind $mats_out docker://xinglab/rmats:v4.1.2 python /rmats/rmats.py";
    $cmd .= " --gtf " . $args->{gtf_file} . " $stranded";
    $cmd .= " --readLength $read_len ";
    $cmd .= " --nthread 4";
    $cmd .= " --tmp $mats_tmp";

    foreach my $c1 ( keys %{$samp->{comps}}) {
        foreach my $c2 ( @{$samp->{comps}->{$c1}} ) {
            print "Running mats for comparison $c1\_vs_$c2\n";
            
            my @rmats = glob("$mats_tmp/*.rmats");
            if(@rmats) {
                unlink @rmats;
            }

            my $result_dir = "$mats_out/$c1\_vs_$c2";
            mkdir $result_dir;

            my $b1_file = "$result_dir/b1.txt";
            my $b2_file = "$result_dir/b2.txt";

            my $b1_content = join(",", find_bams($args->{alignments_dir}, $samp->{group}->{$c1}));
            open(BC1, '>', $b1_file) or die $!;
            print BC1 $b1_content;
            close(BC1);

            my $b2_content = join(",", find_bams($args->{alignments_dir}, $samp->{group}->{$c2}));
            open(BC2, '>', $b2_file) or die $!;
            print BC2 $b2_content;
            close(BC2);

            my $s_cmd = $cmd . " --b1 $b1_file --b2 $b2_file --od $result_dir";
            $s_cmd .= " > $log_dir/mats_$c1\_vs_$c2.log 2>&1";

            print "$s_cmd\n\n";
            system($s_cmd);
        }
    }
}

sub find_bams {
    my $align_dir = shift;
    my $sample_ary = shift;

    my @bams;
    
    foreach my $samp (@{$sample_ary}){
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

    my $mats_out = $args->{output_dir} . "/mats";

    my $cmd = "singularity exec --bind " . $args->{output_dir};
    $cmd .= " --bind " . dirname $args->{ref_fasta};
    $cmd .= " --bind " . dirname $args->{annotation_file};
    $cmd .= " docker://kk988/splicetools:v1.1_1 perl ";

    foreach my $c1 (keys %{$samp->{comps}}) {
        foreach my $c2 (@{$samp->{comps}->{$c1}}) {
            my $mats_dir = "$mats_out/$c1\_vs_$c2";
            my $expr_file = combine_expression_files($samp, $c1, $c2, $mats_dir);
            my $ri_file = "$mats_dir/RI.MATS.JCEC.txt";
            my $se_file = "$mats_dir/SE.MATS.JCEC.txt";
            my $ri_logfile = $args->{output_dir} . "/logs/RIMedley_$c1\_v1_$c2.log";
            my $se_logfile = $args->{output_dir} . "/logs/SEMedley_$c1\_v1_$c2.log";
            my $num_samples = join(",", scalar(@{$samp->{group}->{$c1}}), scalar(@{$samp->{group}->{$c2}}));

            my $ri_medley = $cmd 
                . " /usr/bin/SpliceTools-1.1/bin/RIMedley.pl -r $ri_file"
                . " -a " . $args->{annotation_file}
                . " -g " . $args->{ref_fasta}
                . " -e $expr_file"
                . " -TPM 2,2"
                . " -SN $num_samples"
                . " -f 0.01";
            $ri_medley .= " > $ri_logfile 2>&1";

            print "$ri_medley\n";
            system($ri_medley);

            my $se_medley = $cmd 
                . " /usr/bin/SpliceTools-1.1/bin/SEMedley.pl -s $se_file"
                . " -a " . $args->{annotation_file}
                . " -g " . $args->{ref_fasta}
                . " -e $expr_file"
                . " -TPM 2,2"
                . " -SN $num_samples"
                . " -f 0.01";
            $se_medley .= " > $se_logfile 2>&1";

            print "$se_medley\n\n";
            system($se_medley);

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
    my $final_expr = "$out_dir/expression_file.tsv";
    my $header_line = "target_id";

    my $first_samp = $samp->{samp_list}->[0];
    # First grab column
    my $cmd = "cut -f 1 " . $samp->{ $first_samp }->{kallisto_dir} . "/abundance.tsv > $expr";
    system($cmd);

    # Then iteratively grab the column from each sample's abundance file
    foreach my $sname (@{$samp->{group}->{$grp1}}){
        my $abundance = $samp->{ $sname }->{kallisto_dir} . "/abundance.tsv";
        $header_line .= " $sname";
        my $cmd2 = "awk '{print \$NF}' $abundance | paste $expr - > $tmp_expr; mv $tmp_expr $expr;";
        system($cmd2);
    }

    foreach my $sname (@{$samp->{group}->{$grp2}}){
        my $abundance = $samp->{ $sname }->{kallisto_dir} . "/abundance.tsv";
        $header_line .= " $sname";
        my $cmd2 = "awk '{print \$NF}' $abundance | paste $expr - > $tmp_expr; mv $tmp_expr $expr;";
        system($cmd2);
    }

    # Finally - add the header line:
    my $cmd3 = "echo \"$header_line\" | tr ' ' '\t' > $final_expr; tail -n +2 $expr >> $final_expr;";
    system($cmd3);

    return $final_expr;
}

1;
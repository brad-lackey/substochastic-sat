#!/usr/bin/perl

$queue_size = 2;

(scalar(@ARGV) == 5) or (scalar(@ARGV) == 6) or (scalar(@ARGV) == 8) or die "Usage: ./testrun.pl command <LUT.txt> <filelist.dat> trials tag [seed [\"step weight\" \"runtime\"]]\n";

# Just to be completely clear:
$command = $ARGV[0];
$lut = $ARGV[1];
$filelist = $ARGV[2];
$trials = $ARGV[3];
$tag = $ARGV[4];
$seed = $ARGV[5];
$weight = $ARGV[6];
$runtime = $ARGV[7];

if($seed eq "")
{
    $seed = int(rand(1000000000))
}

printf "SEED: %d..%d, QUEUE_SIZE: %d\n", $seed, $seed+$trials, $queue_size;

# Load int the data file and create the list of all jobs.
open(DAT, "<$filelist");
while(<DAT>){
    chomp;
    my @a = split;
    my $file = $a[0];
    my $optimal = $a[3];
    my $besttime = $a[6];
    for (my $t = $seed; $t < $seed+$trials; $t += 1){
        push(@stack, $command . " $lut $file $optimal $t $weight $runtime");
    }
    $opt{$a[0]} = $optimal;
    $time{$a[0]} = $besttime;
}
close(DAT);

# Randomize the job order with Fisher-Yates.
for (my $i = scalar(@stack)-1; $i>0; $i -= 1) {
    my $r = int rand (1 + $i);
    @stack[$i, $r] = @stack[$r, $i] unless $r == $i;
}

open(LOG, ">$tag.log") || die "Could not open log file.\n";
unlink "$tag.out";

# Run all the jobs.
while ( scalar(@stack) > 0 ){
    $job = shift @stack;
    
    if ( scalar(keys %queue) == $queue_size ){
        my $j = wait();
        print LOG $queue{$j}, "\n";
        delete $queue{$j};
    }
    
    if ( $pid = fork ){
        $queue{$pid} = $job;
    } else {
        @out = `$job`;
        
        $opt = -1;
        $time = -1;
        $loops = -1;
        $updates = -1;
        for $line (@out){
            if ($line =~ /^o/){
                @b = split (/\s/,$line);
                $opt = $b[1];
            }
            if ($line =~ /Walltime/){
                @c = split (/\s/,$line);
                $time = $c[2];
                $loops = $c[4];
                $updates = $c[6];
                # print "$c[1] $c[2] $c[3] $c[4] $c[5] $c[6]\n"
            }
        }
        
        my @a = split(/\s+/,$job);
        open(OUT, ">>$tag.out") || die;
        flock(OUT, 2) || die;
        print OUT "$a[-3] $opt $time $loops $updates\n";
        close(OUT);
        exit(0);
    }
}

# Wait until the last jobs finish.
while( scalar(keys %queue) > 0 ){
    my $j = wait();
    print LOG $queue{$j}, "\n";
    delete $queue{$j};
}
close(LOG);

# Load up the output file and get the results.
open(IN,"<$tag.out");
while( <IN> ){
    chomp;
    my @a = split;
    $count{$a[0]} += 1;
    if ( $a[1] <= $opt{$a[0]} ){
        $hit{$a[0]} += 1;
        if ( $a[2] < $time{$a[0]} ){
            $score{$a[0]} += 1;
        }
        $runtime{$a[0]} += $a[2];
        $loops{$a[0]} += $a[3];
        $updates{$a[0]} += $a[4];
    }

    if ($a[1] == -1){
        $factor{$a[0]} += 1000
    }
    else{
        $factor{$a[0]} += (1+$a[1])/(1+$opt{$a[0]});
    }

    # REMOVE
    # print "$a[0] Factor: $factor{$a[0]}\n";
}
close(IN);

# Produce the final report.
open(REPORT,">$tag.txt");
$count = 0;
$hit = 0;
$score = 0;
$runtime = 0.0;
$loops = 0.0;
$updates = 0;
$factor = 0.0;
for $file (sort keys %count){
    $count += $count{$file};
    print REPORT "$file ";
    $hit += $hit{$file};
    printf REPORT "%d/%d(%.0f%%) ",  $hit{$file}, $count{$file}, (100.0*$hit{$file})/$count{$file};
    $score += $score{$file};
    $runtime += $runtime{$file};
    $loops += $loops{$file};
    $updates += $updates{$file};
    $factor += $factor{$file};
    if ($hit{$file} > 0) {
        printf REPORT "(%d/%d(%.0f%%)) ", $score{$file}, $count{$file}, (100.0*$score{$file})/$hit{$file};
        printf REPORT "%.3fs loops=%.1f updates=%.3f factor=%.3f\n", $runtime{$file}/$hit{$file}, (1.0*$loops{$file})/$hit{$file}, $updates{$file}/$hit{$file}, $factor{$file};
    } else {
        printf REPORT "(%d/%d(%.0f%%)) ", $score{$file}, $count{$file}, 0.0;
        printf REPORT "%.3fs loops=%.1f updates=%.1f factor=%.3f\n", 0.0, 0.0, 0.0, $factor{$file};
    }
}

printf REPORT "Overall: %d/%d(%.0f%%) ", $hit, $count,  (100.0*$hit)/$count;
if ($hit > 0){
    printf REPORT "(%d/%d(%.0f%%)) ", $score, $count, (100.0*$score)/$hit;
    printf REPORT "%.3fs ", $runtime/$hit;
    printf REPORT "%.3f loops ", $loops/$hit;
    printf REPORT "%.3f updates ", $updates/$hit;
    printf REPORT "%.3f factor\n", $factor;
} else {
    printf REPORT "(%d/%d(%.0f%%)) ", $score, $count, 0.0;
    printf REPORT "%.3fs ", 0.0;
    printf REPORT "%.3f loops ", 0.0;
    printf REPORT "%.3f updates ", 0;
    printf REPORT "%.3f factor\n", $factor;
}
close(REPORT);

__END__

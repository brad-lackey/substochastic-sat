#!/usr/bin/perl

$queue_size = 10;

(scalar(@ARGV) == 4) or (scalar(@ARGV) == 6) or die "Usage: ./testrun.pl command <filelist.dat> trials tag [\"step weight\" \"runtime\"]\n";

# Just to be completely clear:
$command = $ARGV[0];
$filelist = $ARGV[1];
$trials = $ARGV[2];
$tag = $ARGV[3];
$weight = $ARGV[4];
$runtime = $ARGV[5];

# Load int the data file and create the list of all jobs.
open(DAT, "<$filelist");
while(<DAT>){
  chomp;
  my @a = split;
  my $file = $a[0];
  my $optimal = $a[3];
  my $besttime = $a[6];
  for (my $t = 0; $t < $trials; $t += 1){
    push(@stack, $command . " $file $optimal $t $weight $runtime");
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
    for $line (@out){
      if ($line =~ /^o/){
	@b = split (/\s/,$line);
	$opt = $b[1];
      }
      if ($line =~ /Walltime/){
	@c = split (/\s/,$line);
	$time = $c[2];
	$loops = $c[4];
      }
    }

    my @a = split(/\s+/,$job);
    open(OUT, ">>$tag.out") || die;
    flock(OUT, 2) || die;
    print OUT "$a[1] $opt $time $loops\n";
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
  $runtime{$a[0]} += $a[2];
  $loops{$a[0]} += $a[3];
  if ( $a[1] == $opt{$a[0]} ){
    $hit{$a[0]} += 1;
    if ( $a[2] < $time{$a[0]} ){
      $score{$a[0]} += 1;
    }
  }
}
close(IN);

# Produce the final report.
open(REPORT,">$tag.txt");
$count = 0;
$hit = 0;
$score = 0;
$runtime = 0.0;
for $file (sort keys %count){
	$count += $count{$file};
  print REPORT "$file ";
	$hit += $hit{$file};
  printf REPORT "%d/%d(%.0f%%) ",  $hit{$file}, $count{$file}, (100.0*$hit{$file})/$count{$file};
	$score += $score{$file};
  printf REPORT "(%d/%d(%.0f%%)) ", $score{$file}, $count{$file}, (100.0*$score{$file})/$count{$file};
	$runtime += $runtime{$file};
  printf REPORT "%.3fs loops=%.1f\n", $runtime{$file}/$count{$file}, $loops{$file}/$count{$file};
}
printf REPORT "Overall: %d/%d(%.0f%%) ", $hit, $count,  (100.0*$hit)/$count;
printf REPORT "(%d/%d(%.0f%%)) ", $score, $count, (100.0*$score)/$count;
printf REPORT "%.3fs\n", $runtime/$count;
close(REPORT);

__END__

#!/usr/bin/perl

$queue_size = 7;

(scalar(@ARGV) == 3) or (scalar(@ARGV) == 5) or die "Usage: ./testrun.pl <data-file> trials tag [\"step weight\" \"runtime\"]\n";

open($dat, "<", $ARGV[0]);

# Load int the data file and create the list of all jobs.
while(<$dat>){
  chomp;
  @a = split;
  for ($t = 0; $t < $ARGV[1]; $t += 1){
    push(@stack, "./ssmc $a[0] $a[3] $t $ARGV[3] $ARGV[4]");  #<--------- Here is where the list of command is created.
  }
  $opt{$a[0]} = $a[3];
  $time{$a[0]} = $a[6];
}

# Randomize the job order with Fisher-Yates.
for (my $i = scalar(@stack)-1; $i>0; $i -= 1) {
  my $r = int rand (1 + $i);
  @stack[$i, $r] = @stack[$r, $i] unless $r == $i;
}

open(LOG, ">$ARGV[2].log") || die "Could not open log file.\n";
unlink "$ARGV[2].out";

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

    @a = split(/\s+/,$job);
    open(OUT, ">>$ARGV[2].out") || die;
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

open(IN,"<$ARGV[2].out");

# Load up the output file and get the results.
while( <IN> ){
  chomp;
  @a = split;
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

open(REPORT,">$ARGV[2].txt");

# Produce the final report.
for $file (sort keys %count){
  print REPORT "$file ";
  printf REPORT "%d/%d(%.0f%%) ",  $hit{$file}, $count{$file}, (100.0*$hit{$file})/$count{$file};
  printf REPORT "(%d/%d(%.0f%%)) ", $score{$file}, $count{$file}, (100.0*$score{$file})/$count{$file};
  printf REPORT "%.3fs loops=%.1f\n", $runtime{$file}/$count{$file}, $loops{$file}/$count{$file};
}

close(REPORT);

__END__

#!/usr/bin/perl

open(SAT,"<$ARGV[0]");
while(<SAT>){
  if ( /^p/ ){
    @pline = split;
    $vars = int($pline[2]);
    $clss = int($pline[3]);
    if ( $pline[1] =~ /wcnf/ ){
      $mode = 0;
    }
    else {
      $mode = 1;
    }
    last;
  }
}

push @list, 0;
while(<SAT>){
  chomp;
  @cline = split(/\s+/);
  if ( $mode == 1 ){
    push @list, 1;
  }
  push @list, @cline;
  
}
close(SAT);

#print join (" ", @list);
#print "\n";

open(OUT,"<$ARGV[1]");
while(<OUT>){
  if ( /^o/ ){
    print;
  }
  if ( /^v/ ) {
    chomp;
    @v = split(/\s+/);
    $o = 0;
    $wt = 0;
    $valid = 1;
    for ($j=0; $j<scalar(@list); $j+=1){
      if ($list[$j] == 0){
	if ($valid == 0 ){
	  $o += $wt;
	}
	$valid = 0;
	$j += 1;
	$wt = $list[$j];
      }
      else {
	if ( $list[$j] < 0 ) {
	  if ( $v[-$list[$j]] < 0 ){
	    $valid = 1;
	  }
	}
	else {
	  if ( $v[$list[$j]] > 0 ){
	    $valid = 1;
	  }	  
	}
      }
    }
    print "v $o\n";
  }
}
close(OUT);

__END__

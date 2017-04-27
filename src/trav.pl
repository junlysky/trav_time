#!/usr/bin/perl
use strict;
# perl code for computing first arrival times in layered velocity model.
# by junlysky.info@gmail.com

my $r0 = 6371.;         # for flap convert 
my $r_depth = 0.;	# station receiver depth.
my $deg2km=1;
my $flat=0;		    # Earth flattening transformation.
my $kappa = 0;		# the input model 3rd column is vp, not vp/vs ratio.
my $rdep = "";

###################################### end of defaut parameters

@ARGV >= 2 or die "Usage:
      trav.pl -Mmodel/depth[/f_or_k] [-Rrdep] distance
      -R: the station receiver depth ($r_depth=0)
 \n";


my ($model, $s_depth);


foreach (grep(/^-/,@ARGV)) {
   my $opt = substr($_,1,1);
   my @value = split(/\//,substr($_,2));
   if ($opt eq "D") {
    	$deg2km = 6371*3.14159/180.;
   	} 
   elsif ($opt eq "M") {
   		$model = $value[0];
    	$s_depth = $value[1] if $#value > 0;
        $flat = 1 if $value[2] eq "f";
        $kappa = 1 if $value[2] eq "k";
   }   
   elsif ($opt eq "R") {
     $r_depth = $value[0];
     $rdep = "_$r_depth";
   } 

   else {
     print STDERR "Error **** Wrong options\n";
     exit(0);
   }
}
my (@dist) = grep(!/^-/,@ARGV);


if ($flat) {
   $s_depth = $r0*log($r0/($r0-$s_depth));
   $r_depth = $r0*log($r0/($r0-$r_depth));
}

# input velocity model
my (@th, @vs, @vp, @rh, @qa, @qb);
my ($src_layer, $rcv_layer);
&read_model();
print " thickness: @th \n vs: @vs\n vp: @vp\n rh: @rh\n qa: @qa\n qb: @qb\n";
my $freeSurf = $th[0]>0. || $#th<1 ? 1 : 0;
if ( $freeSurf && ($s_depth<0. || $r_depth<0.) ) {
   print STDERR "Error **** The source or receivers are located in the air!!!???\n";
   exit(0);
}
if ($s_depth<$r_depth) {
   $src_layer = insert_intf($s_depth);
   $rcv_layer = insert_intf($r_depth);
} else {
   $rcv_layer = insert_intf($r_depth);
   $src_layer = insert_intf($s_depth);
}

if ( $vp[$src_layer] != $vp[$src_layer-1] ) {
   print STDERR "Error **** The source is located at a real interface\n";
   exit(0);
}
my $num_layer = $#th + 1;

print " station layer is: $rcv_layer\n source layer is: $src_layer \n";
# compute first arrival times
my (%t0, %sac_com);
&ps_arr();

sub read_model {
   open(MODEL,"$model") or die "could not open $model\n";
   my $fl = 1.;
   my $r = $r0;
   my $i = 0;
   while (<MODEL>){
      ($th[$i], $vs[$i], $vp[$i], $rh[$i], $qb[$i], $qa[$i]) = split;
      $r -= $th[$i];
      $fl = $r0/($r + 0.5*$th[$i]) if $flat;
      $th[$i] *= $fl;
      $vs[$i] *= $fl;
      if ($kappa) {
         $vp[$i] = $vp[$i]*$vs[$i];	# 3rd column is Vp/Vs
      } else {
         $vp[$i] *= $fl;
      }
      if (!$rh[$i] or $rh[$i] > 20.) {	# 4th column is Qs
         $qb[$i] = $rh[$i];
         $rh[$i] = 0.77 + 0.32*$vp[$i];
      }
      $qb[$i] = 500 unless $qb[$i];
      $qa[$i] = 2*$qb[$i] unless $qa[$i];
      $i++;
   }
   close(MODEL);
   $th[$i-1] = 0.;
}

# compute the source layer
sub insert_intf {
   my ($zs) = @_;
   my $n = $#th + 1;
   my $dep = 0.;
   my $i = 0;
   for(; $i<$n; $i++) {
      $dep += $th[$i];
      last if $dep>$zs || $i==$n-1;
   }
   my $intf = $i;
   return $intf if ($i>0 && $zs==$dep-$th[$i]) || ($i==0 && $zs==0.);
   my $dd = $dep-$zs;
   for($i=$n; $i>$intf; $i--) {
       $th[$i] = $th[$i-1];
       $vs[$i] = $vs[$i-1];
       $vp[$i] = $vp[$i-1];
       $rh[$i] = $rh[$i-1];
       $qa[$i] = $qa[$i-1];
       $qb[$i] = $qb[$i-1];
   }
   $th[$intf]  -= $dd;
   $th[$intf+1] = $dd if $dd>0.;
   if ($th[0]<0.) {
      $s_depth -= $th[0];
      $r_depth -= $th[0];
      $th[0]=0.;
   }
   return $intf+1;
}

sub ps_arr {
   my ($j, $tp, $ts, $pa, $sa, $dn, @aaa);

   # calculate arrival time for P and S
   open(TRAV2,"| trav2 > junk.p") or die "couldn't run trav2\n";
 #  print "$num_layer,$src_layer,$rcv_layer \n";

   printf TRAV2 "%d %d %d\n", $num_layer,$src_layer,$rcv_layer;
   for ($j=0;$j<$num_layer;$j++) {
      printf TRAV2 "%11.4f %11.4f\n",$th[$j],$vp[$j];
   }
   printf TRAV2 "%d\n",$#dist + 1;
   foreach (@dist) {
       printf TRAV2 "%10.4f\n",$_*$deg2km;
   }
   close(TRAV2);
   open(TRAV2,"| trav2 > junk.s") or die "couldn't run trav2\n";
   printf TRAV2 "%d %d %d\n", $num_layer,$src_layer,$rcv_layer;
   for ($j=0;$j<$num_layer;$j++) {
      printf TRAV2 "%11.4f %11.4f\n",$th[$j],$vs[$j];
   }
   printf TRAV2 "%d\n",$#dist + 1;
   foreach (@dist) {
       printf TRAV2 "%10.4f\n",$_*$deg2km;
   }
   close(TRAV2);

   open(TRAV2,"paste junk.p junk.s |");
   my $vps="vps.dat";
   open(VPS," > $vps");
   $j=0;
   while (<TRAV2>) {
     @aaa = split;
     $tp = $aaa[1];
     if ($aaa[2]>$tp && $aaa[3]<1/7.) { # down going Pn
        $pa = $vp[$src_layer]*$aaa[3]; $dn=1;
     } else {
        $pa = $vp[$src_layer]*$aaa[4]; $dn=-1;
     }
     $pa = atan2($pa,$dn*sqrt(abs(1.-$pa*$pa)))*180/3.14159;
     $ts = $aaa[6];
     if ($aaa[7]>$ts && $aaa[8]<1/4.) { # down going Sn
        $sa = $vs[$src_layer]*$aaa[8]; $dn=1;
     } else {
        $sa = $vs[$src_layer]*$aaa[9]; $dn=-1;
     }
     $sa = atan2($sa,$dn*sqrt(abs(1.-$sa*$sa)))*180/3.14159;
#     $sac_com{$dist[$j]}=sprintf("t1 $tp t2 $ts user1 %6.2f user2 %6.2f",$pa,$sa);
     print VPS "$dist[$j] $s_depth     tp $tp    ts $ts\n";
     $t0{$dist[$j]} = $tp;
#     print STDERR "$sac_com{$dist[$j]} $t0{$dist[$j]}\n";
      print STDERR "$dist[$j] $s_depth   tp $tp    ts $ts\n";
     $j++;
   }
   close (TRAV2);
   close (VPS);
}

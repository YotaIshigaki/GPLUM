#!/usr/bin/perl -w
use strict;

# This script needs one argument from the command line.

# Check the number of the arguments
if ( @ARGV != 1 ) {
   die "The number of the arguments should be one.\n";
}
# Check the value of the second argument.
my $restart = $ARGV[0];
my $mktable;
if ( "$restart" eq "restart_run" ) {
   $restart = 1;
   $mktable = 0;
} elsif ( "$restart" eq "initial_run" ) {
   $restart = 0;
   $mktable = 0;
} elsif ( "$restart" eq "make_tables" ) {
   $restart = 0;
   $mktable = 1;
} else {
   die "The content of the second argument is invalid\n";
}

#-----------------------------------------------------------
# [1] Declare preprocess keywords and those values by hash
#-----------------------------------------------------------
our $pointer=-1;
our %pp_keywords;
our @ordering;

our $pointer_wc=-1;
our @pp_keywords_wc;
our @pp_values_wc;
our @pp_conditions_wc;
# `wc` is a short of `with conditions`

our $pointer_err=-1;
our @err_list;

#-----------------------------------------------------------
# [2] Define keywords and values
#-----------------------------------------------------------
# [!!!!! Important Notes !!!!!]
#   - Never change the keywords and values marked by #!.
#   - Macros with `future` are reserved macros for future implementation.
#

#===========================================
#  Base settings
#===========================================
&set_pp_keywords("Restart_Mode"   ,"$restart"); #!
&set_pp_keywords("Make_Rate_Tables" ,"$mktable"); #!
&set_pp_keywords("Simulation_Dimension","3"); #!

# PM Test options
&set_pp_keywords("NaCl_Crystal_Mode","1");
&set_pp_keywords("Charge_Neutrality","0");
&set_pp_keywords("Self_Interaction_Term","1");

#-----------------------------------------------------------
# [3] Open file & Output
#-----------------------------------------------------------
open File, ">preprocess_keywords.h" or die "cannot make the file.";
   # [3-1] write simple macro definitions
   if ($pointer ne "-1") {
      for (my $i=0;$i<=$pointer;$i++) {
         my $output = "#define " . $ordering[$i] . " (" . $pp_keywords{$ordering[$i]} . ")";
         print File "$output\n";
      }
   }
   print File "\n\n";
   # [3-2] write #pragma once
   print File "#pragma once\n";
   # [3-3] write conditionally-defined macros
   print File "//=============================\n";
   print File "// Conditionally-defined macros\n";
   print File "//=============================\n";
   if ($pointer_wc ne "-1") {
      for (my $i=0;$i<=$pointer_wc;$i++) {
         my $output = "#if (" . $pp_conditions_wc[$i] . ")\n"
	 	    . "#ifdef " . $pp_keywords_wc[$i] . "\n"
		    . "#undef " . $pp_keywords_wc[$i] . "\n"
		    . "#endif\n"
                    . "#define " . $pp_keywords_wc[$i] . " (" . $pp_values_wc[$i] . ")\n"
                    . "#endif";
         print File "$output\n\n";
      }
   }
   print File "\n\n";
   # [3-4] write error detect script
   print File "//=============================\n";
   print File "// Error conditions            \n";
   print File "//=============================\n";
   if ($pointer_err ne "-1") {
      for (my $i=0;$i<=$pointer_err;$i++) {
         my $output = "#if (" . $err_list[$i] . ")\n"
                    . "#error\n"
                    . "#endif";
         print File "$output\n\n";
      }
   }
close File;


sub set_pp_keywords (\$\$) {
  my $key   = shift;
  my $value = shift;

  $pointer++;
  $pp_keywords{$key} = $value;
  $ordering[$pointer] = $key;

}

# Conditionally overwrite
sub conditionally_overwrite (\$\$\$) {
  my $key       = shift;
  my $value     = shift;
  my $condition = shift;

  $pointer_wc++;
  $pp_keywords_wc[$pointer_wc]   = $key;
  $pp_values_wc[$pointer_wc]     = $value;
  $pp_conditions_wc[$pointer_wc] = $condition;

}

# Error Check
sub error_condition (\$) {
  my $err_condition = shift;

  $pointer_err++;
  $err_list[$pointer_err] = $err_condition;

}

#!/usr/bin/perl -w
use strict;

# (1) Fortran files
my @target_files = glob "*.h *.cpp *.hpp";

# (2) Macros that defined in create_preprocess_keywords.pl
our %defined_macros;
open File, "<create_preprocess_keywords.pl" or die "file is not found.";
while(defined (my $line = <File>)){
   chomp $line;
   if ($line =~ /&set_pp_keywords\("(\w+)"/) {
      #print "$1\n";
      $defined_macros{$1}=1;
   }
}
close(File);


# (3) Check undefined macros are used or not
while(defined (my $file_name = pop(@target_files))){
   print "Checking --- $file_name\n";

   # (3-1) Detect preprocessor instruction
   my @used_macros = ();
   open File, "<$file_name" or die "specified source file is not found.";
      while(defined (my $line = <File>)){
	      if (my @tmp = &get_used_macros($line)) {
	         push (@used_macros, @tmp);
            #print "@tmp\n";
 	      }
      }
   close(File);

   # (3-2) Choose Unique macros
   my %tmp;
   my @unique_macros = sort (grep(  !$tmp{$_}++, @used_macros ));

	# (3-2) Check the validity of the unique macro
   if ( &check_validity(@unique_macros) eq 1){
      exit(1);
   }
}

#----------------------------------------------------------------------
# SUBROUTINE <output_list>
#----------------------------------------------------------------------
sub output_list (\@) {

   while (my $element = pop(@_)) {
      print "$element\n";
   }
}

#----------------------------------------------------------------------
# SUBROUTINE <check_validity>
#----------------------------------------------------------------------
sub check_validity (\@) {
   my $err_flag = 0;

   while (defined (my $macro = pop(@_))) {
      if (!defined($defined_macros{$macro})) { 
         print "$macro is not defined in create_preprocess_keywords.pl\n";
         $err_flag = 1;
      }  
   }

   if ($err_flag eq 0) {
      return 0;
   } else {
      return 1;
   }

}

#----------------------------------------------------------------------
# SUBROUTINE <get_used_macros>
#----------------------------------------------------------------------
sub get_used_macros(\$) {
   my $line = shift;
   my @macros_in_line;

   if ($line =~ /(?:#if|#elif)\s*\((.*)\)/ ){
      if (!($` =~ /!/)) { # skip, if Fortran comment character is containd.
         @macros_in_line = $1 =~ /(\w{2,})/g;
      }
   }

   return @macros_in_line;
}

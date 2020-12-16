#!/usr/bin/ruby

=begin
This program creates a log file of makeflags.
=end

$stderr.print "Make RunLogsAG.h, Start!\n"

filename = ARGV[0]

read_file = open(filename)
output_file = open("RunLogsAG.h","w")

output_file.print "#include \"config.h\"\n"
output_file.print 'static void WriteMakeflagsAG(char Fname[]){',"\n"
output_file.print "\n"
output_file.print "    FILE *fp;\n"
output_file.print "    FileOpen(fp,Fname,\"w\");\n"
output_file.print "\n"
output_file.print '    fprintf(fp,"# Auto generation in %s\\n",__DATE__);',"\n"
output_file.print "\n"
    while text = read_file.gets do
        text.chomp!
        text.gsub!('"','\"')
        output_file.print '    fprintf(fp,"'+text+'\\n");',"\n"
#print '    fprintf(stderr,"'+text+'\\n");',"\n"
    end

output_file.print "\n"
output_file.print "    fclose(fp);\n"
output_file.print "\n"
output_file.print '    return ;',"\n"
output_file.print '}',"\n"

read_file.close 
output_file.close 

$stderr.print "Make RunLogsAG.h, Finished!\n"

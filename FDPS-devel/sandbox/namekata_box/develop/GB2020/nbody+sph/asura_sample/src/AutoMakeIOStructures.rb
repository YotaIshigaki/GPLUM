#!/usr/bin/ruby

=begin
This program analyzes DataStructure.h and creates functions for data output.
Structures of Pbody, Phydro and Pstar are analyzed.
=end

filename = "./DataStructures.h"

##
## time stamp をみて、ファイルが新しいようだったら更新する。 
##

def return_array_size(str)
    if /\[/ =~ str
        /\]/ =~ $'
        array_size = $`.to_i
        return array_size
    else
        return 1
    end
end

def return_array_name(str)
    if /\[/ =~ str
        return $`
    else
        /;/ =~ str
        return $`
    end
end

def return_type_name_size(str)
    str.chomp!
    str.strip!
    column = str.split(/\s+/)

    if /unsigned/ =~ column[0]
        if /long/ =~ column[1]
            if /int/ =~ column[2]
                array_size = return_array_size(column[3])
                array_name = return_array_name(column[3])    
                return column[0]+' '+column[1]+' '+column[2],array_name,array_size
            else
                array_size = return_array_size(column[2])
                array_name = return_array_name(column[2])
                return column[0]+' '+column[1],array_name,array_size
            end
        elsif /int/ =~ column[1]
            array_size = return_array_size(column[2])
            array_name = return_array_name(column[2])
            return column[0]+' '+column[1],array_name,array_size
        end
    else
        array_size = return_array_size(column[1])
        array_name = return_array_name(column[1])
        return column[0],array_name,array_size
    end
end

def write_contains_of_structure(current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')
        subtext.gsub!('*','\*')
        if /<S,C>/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',type,' ',array_name,';'"\n"
            else
                current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
            end
        elsif /<S>/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',type,' ',array_name,';'"\n"
            else
                current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
            end
        elsif /^};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_compact(current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')
        subtext.gsub!('*','\*')
        if /<S,C>/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',type,' ',array_name,';'"\n"
            else
                current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
            end
        elsif /^};/ =~ subtext
            break
        end
    end
end

def write_definition_of_structure(structure_name,current_read_file,current_output_file,type)

        current_output_file.print structure_name,"\n"
        current_output_file.print "\n"

        if type === "standard"
            write_contains_of_structure(current_read_file,current_output_file)
        elsif type === "compact"
            write_contains_of_structure_for_compact(current_read_file,current_output_file)
        else
            STDOUT.print "File type error!\n"
        end

        current_output_file.print "\n"
        current_output_file.print '};',"\n"
        current_output_file.print "\n"

end


def write_contains_of_structure_for_copy(structure_name,source_structure_name,current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')
        subtext.gsub!('*','\*')
        if /<S,C>/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    source_structure_name,'[index]->',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /<S>/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    source_structure_name,'[index]->',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^};/ =~ subtext
            break
        end
    end
end

def read_contains_of_structure_from_copy(source_structure_name,dest_structure_name,current_read_file,current_output_file,version)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')
        subtext.gsub!('*','\*')
        if /<S,C>/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',dest_structure_name,'.',array_name,' = ',
                    source_structure_name,'.',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',dest_structure_name,'.',array_name,'[',i,']'' = ',
                        source_structure_name,'.',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /<S>/ =~ subtext
            if version === "standard" then
                type,array_name,array_size = return_type_name_size(subtext)
                if array_size == 1
                    current_output_file.print '    ',dest_structure_name,'.',array_name,' = ',
                        source_structure_name,'.',array_name,';'"\n"
                else
                    array_size.times {|i|
                        current_output_file.print '    ',dest_structure_name,'.',array_name,'[',i,']'' = ',
                            source_structure_name,'.',array_name,'[',i,']'';'"\n"
                    }
                end
            end
        elsif /^};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_copy_compact(structure_name,source_structure_name,current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')
        subtext.gsub!('*','\*')
        if /<S,C>/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    source_structure_name,'[index]->',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^};/ =~ subtext
            break
        end
    end
end

def read_contains_of_structure_from_copy_compact(structure_name,source_structure_name,current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')
        subtext.gsub!('*','\*')
        if /<S,C>/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',source_structure_name,'.',array_name,' = ',
                    structure_name,'.',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',source_structure_name,'.',array_name,'[',i,']'' = ',
                        structure_name,'.',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^};/ =~ subtext
            break
        end
    end
end

def write_function_for_copy(function_name,source_structure_name,current_read_file,current_output_file,type)

        current_output_file.print function_name,"\n"
        current_output_file.print "\n"

        split_structure_name = function_name.split(/\s+/)
        current_structure_name = split_structure_name[1]

        current_output_file.print "\t"+split_structure_name[0]+' '+split_structure_name[1]+";\n"
        current_output_file.print "\n"

        if type === "standard"
            write_contains_of_structure_for_copy(current_structure_name,
                    source_structure_name,current_read_file,current_output_file)
        elsif type === "compact"
            write_contains_of_structure_for_copy_compact(current_structure_name,
                    source_structure_name,current_read_file,current_output_file)
        else
            STDOUT.print "File type error!\n"
        end

        current_output_file.print "\n"
        current_output_file.print "\treturn "+split_structure_name[1]+";\n"

        current_output_file.print "\n"
        current_output_file.print '};',"\n"
        current_output_file.print "\n"

end

def read_function_from_copy(structure_name,current_read_file,current_output_file,version)

        type_of_structure = "Struct"+structure_name
        dest_structure_name = structure_name+"Temporary"
        source_structure_name = structure_name+"IO"
        
        if version === "standard"
            function_name = type_of_structure+" "+"CopyTemporalStructureTo"+
                structure_name+"("+type_of_structure+"IO"+" "+source_structure_name+"){"
        elsif version === "compact"
            function_name = type_of_structure+" "+"CopyTemporalStructureTo"+
                structure_name+"Compact"+"("+type_of_structure+"IO"+"Compact"+" "+source_structure_name+"Compact){"
            source_structure_name = structure_name+"Compact"
        else
            STDOUT.print "File type error!\n"
        end

        current_output_file.print function_name,"\n"
        current_output_file.print "\n"

        current_output_file.print "\t"+type_of_structure+' '+dest_structure_name+";\n"
        current_output_file.print "\n"

        read_contains_of_structure_from_copy(source_structure_name,
            dest_structure_name,current_read_file,current_output_file,version)

        current_output_file.print "\n"
        current_output_file.print "\treturn "+dest_structure_name+";\n"

        current_output_file.print "\n"
        current_output_file.print '};',"\n"
        current_output_file.print "\n"

end

=begin
First, this script writes "Standard Output" members in structures.
Next, the script writes "Compact Output" members in the structures.
At last, generate functions, which is copying from ordinary structure to temporal structures.
=end

$stderr.print "Make StructIOAG.h, Start!\n"
read_file = open(filename)
output_header = open("StructIOAG.h","w")
output_header.print "\n"
output_header.print '// Auto generation in ',Time.now,"\n"
output_header.print "\n"

output_body = open("StructIOAG.c","w")
output_body.print "#define \"config.h\"\n"
output_body.print "#define \"StructIOAG.h\"\n"
output_body.print "\n"
output_body.print '// Auto generation in ',Time.now,"\n"
output_body.print "\n"

while text = read_file.gets do
    text.chomp!
    text.gsub!('"','\"')
    case text
    when "struct StructPbody_tag{"

        filepointer = read_file.pos
        write_definition_of_structure("struct StructPbodyIO{",read_file,output_header,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPbodyIOCompact{",read_file,output_header,"compact")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPbodyIO CopyPbodyToTemporalStructure(const int index){",
                "Pbody",read_file,output_body,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPbodyIOCompact CopyPbodyToTemporalStructureCompact(const int index){",
                "Pbody",read_file,output_body,"compact")


        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pbody",read_file,output_body,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pbody",read_file,output_body,"compact")

    when "struct StructPhydro_tag{"

        filepointer = read_file.pos
        write_definition_of_structure("struct StructPhydroIO{",read_file,output_header,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPhydroIOCompact{",read_file,output_header,"compact")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPhydroIO CopyPhydroToTemporalStructure(const int index){",
                "Phydro",read_file,output_body,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompact(const int index){",
                "Phydro",read_file,output_body,"compact")

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Phydro",read_file,output_body,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Phydro",read_file,output_body,"compact")

    when "struct StructPstar_tag{"

        filepointer = read_file.pos
        write_definition_of_structure("struct StructPstarIO{",read_file,output_header,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPstarIOCompact{",read_file,output_header,"compact")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPstarIO CopyPstarToTemporalStructure(const int index){",
                "Pstar",read_file,output_body,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPstarIOCompact CopyPstarToTemporalStructureCompact(const int index){",
                "Pstar",read_file,output_body,"compact")

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pstar",read_file,output_body,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pstar",read_file,output_body,"compact")

    when "struct StructPsink_tag{"

        filepointer = read_file.pos
        write_definition_of_structure("struct StructPsinkIO{",read_file,output_header,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPsinkIOCompact{",read_file,output_header,"compact")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPsinkIO CopyPsinkToTemporalStructure(const int index){",
                "Psink",read_file,output_body,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompact(const int index){",
                "Psink",read_file,output_body,"compact")

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Psink",read_file,output_body,"standard")
        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Psink",read_file,output_body,"compact")
    end
end

read_file.close 
output_header.close 
output_body.close 

$stderr.print "Make StructIOAG.h, Finished!\n"


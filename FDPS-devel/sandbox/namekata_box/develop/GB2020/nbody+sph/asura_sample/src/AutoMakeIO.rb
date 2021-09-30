#!/usr/bin/ruby

=begin
This program analyzes DataStructure.h and creates functions for data output.
Structures of Pbody, Phydro and Pstar are analyzed.
=end

filename = "./DataStructures.h"

# Check number of elements used in CELib
def return_element_size
    hostname = ENV["HOSTNAME"]
    puts hostname

    flag = 0
    targetMakeFile = ""
    # Check Makefile
    File.open('./machines/Makefile.switch') do |f|
        f.each_line do |line|
            if flag == 1
                temp = line.split(" ")
                targetMakeFile = temp[1]
                #TargetMakeFile = line.split!
                break
            end
            if line.include?(hostname)
                flag = 1
            end
        end
    end
    puts targetMakeFile

    #get uname
    username = ENV['USER']
    puts username

    targetPath = "/home/hiraiyt/work/celib/src"
    File.open(targetMakeFile) do |f|
        f.each_line do |line|
            if /^INCLUDEPATH/ =~ line
            #if line.include?("INCLUDEPATH")
            if line.include?("CELib")
            if line.include?(username)
                temp = line.split("-I")
                targetPath = temp[1]
            end
            end
            end
        end
    end
    targetPath.strip!
    puts targetPath

    targetHeader = targetPath + "/CELib.h"
    puts targetHeader
    fstart = 0
    fend = 0
    File.open(targetHeader) do |f|
        f.each_line do |line|
            puts line
            if line.include?("CELibYield_H,")
                fstart = f.lineno
            end
            if line.include?("CELibYield_Number")
                fend = f.lineno
                break
            end
        end
    end
    puts fend-fstart 
    return fend-fstart
end

Size = return_element_size
puts Size
#exit

##
## time stamp をみて、ファイルが新しいようだったら更新する。 
##

def return_array_size(str)
    str_copy = str

    #  a = str.count("[")
    #  multiarraysize = Array.new(3,0)
    #  if (a>1)
        #  for i in 0..a do
            #  top = str.index("[",i)
            #  bottom = str.index("]",i)
            #  multiarraysize[i] = str[top+1..bottom-1]
            #  #b = str[top+1..bottom-1]
            #  puts multiarraysize[i]
        #  end
        #  #return multiple?, array
        #  return true, multiarraysize
    #  else
        if /\[/ =~ str
            /\]/ =~ $'
            body = $`

            if /CELibYield_Number/ =~ body
                return Size
            end

            array_size = body.to_i
            if array_size.is_a?(Integer)
                #puts str
                #puts body
                return array_size
            end

        else
            return 1
        end
    #  end
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
            if /long/ =~ column[2]
                if /int/ =~ column[3]
                    array_size = return_array_size(column[4])
                    array_name = return_array_name(column[4])    
                    return column[0]+' '+column[1]+' '+column[2]+' '+column[3],array_name,array_size
                else
                    array_size = return_array_size(column[3])
                    array_name = return_array_name(column[3])    
                    return column[0]+' '+column[1]+' '+column[2],array_name,array_size
                end
            elsif /int/ =~ column[2]
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
        # if array_size == Size
            # puts "++++++++++++++++++++++++++++++++++"
            # puts array_size
            # puts array_name
        # end
        return column[0],array_name,array_size
    end
end

def write_contains_of_structure(current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp! #
        subtext.strip! #
        subtext.gsub!('"','\"') #

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(bool|int|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',type,' ',array_name,';'"\n"
            else
                current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
            end
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_compact(current_read_file,current_output_file)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /^\};/ =~ subtext
            break
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        else
            type,array_name,array_size = return_type_name_size(subtext)
=begin
            if array_size == 1
                current_output_file.print '    ',type,' ',array_name,';'"\n"
            else
                current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
            end
=end
            if array_size == 1 && /double/ =~ type
                current_output_file.print '    float ',array_name,';'"\n"
            elsif array_size == 1
                current_output_file.print '    ',type,' ',array_name,';'"\n"
            elsif /double/ =~ type
                current_output_file.print '    float ',array_name,'[',array_size,'];'"\n"
            else
                current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
            end
        end
    end
end

def write_contains_of_structure_for_compact_double(current_read_file,current_output_file)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /^\};/ =~ subtext
            break
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        else
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',type,' ',array_name,';'"\n"
            else
                current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
            end
        end
    end
end

def write_contains_of_structure_for_compact_mix(current_read_file,current_output_file)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /^\};/ =~ subtext
            break
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        else
            if /<MIX>/ =~ subtext
                type,array_name,array_size = return_type_name_size(subtext)
                if array_size == 1
                    current_output_file.print '    ',type,' ',array_name,';'"\n"
                else
                    current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
                end
            else
                type,array_name,array_size = return_type_name_size(subtext)
                if array_size == 1 && /double/ =~ type
                    current_output_file.print '    float ',array_name,';'"\n"
                elsif array_size == 1
                    current_output_file.print '    ',type,' ',array_name,';'"\n"
                elsif /double/ =~ type
                    current_output_file.print '    float ',array_name,'[',array_size,'];'"\n"
                else
                    current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
                end
            end
        end
    end
end

def write_contains_of_structure_for_lean(current_read_file,current_output_file)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            $stderr.print subtext+"\n"
            current_output_file.print subtext,"\n"
        elsif /^\};/ =~ subtext
            break
        elsif /<LEAN>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        else
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1 && /double/ =~ type
                current_output_file.print '    float ',array_name,';'"\n"
            elsif array_size == 1
                current_output_file.print '    ',type,' ',array_name,';'"\n"
            elsif /double/ =~ type
                current_output_file.print '    float ',array_name,'[',array_size,'];'"\n"
            else
                current_output_file.print '    ',type,' ',array_name,'[',array_size,'];'"\n"
            end
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
        elsif type === "compact_double"
            write_contains_of_structure_for_compact_double(current_read_file,current_output_file)
        elsif type === "compact_mix"
            write_contains_of_structure_for_compact_mix(current_read_file,current_output_file)
        elsif type === "lean"
            write_contains_of_structure_for_lean(current_read_file,current_output_file)
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
    end
end

def read_contains_of_structure_from_copy(source_structure_name,dest_structure_name,current_read_file,current_output_file,version)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',dest_structure_name,'.',array_name,' = ',
                    source_structure_name,'.',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',dest_structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'.',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def read_contains_of_structure_from_copy_lean(source_structure_name,dest_structure_name,current_read_file,current_output_file,version)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<LEAN>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',dest_structure_name,'.',array_name,' = ',
                    source_structure_name,'.',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',dest_structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'.',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_copy_compact(structure_name,source_structure_name,current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)

            if array_size == 1 && /double/ =~ type
                current_output_file.print '    ',structure_name,'.',array_name,' = (float)',
                    source_structure_name,'[index]->',array_name,';'"\n"
            elsif array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    source_structure_name,'[index]->',array_name,';'"\n"
            elsif /double/ =~ type
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = (float)',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            end
=begin
            if array_size == 1 
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    source_structure_name,'[index]->',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            end
=end
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_copy_compact_double(structure_name,source_structure_name,current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)

            if array_size == 1 && /double/ =~ type
                current_output_file.print '    ',structure_name,'.',array_name,' = (double)',
                    source_structure_name,'[index]->',array_name,';'"\n"
            elsif array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    source_structure_name,'[index]->',array_name,';'"\n"
            elsif /double/ =~ type
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = (double)',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_copy_compact_mix(structure_name,source_structure_name,current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            if /<MIX>/ =~ subtext
                type,array_name,array_size = return_type_name_size(subtext)
                if array_size == 1 && /double/ =~ type
                    current_output_file.print '    ',structure_name,'.',array_name,' = (float)',
                        source_structure_name,'[index]->',array_name,';'"\n"
                elsif array_size == 1
                    current_output_file.print '    ',structure_name,'.',array_name,' = ',
                        source_structure_name,'[index]->',array_name,';'"\n"
                elsif /double/ =~ type
                    array_size.times {|i|
                        current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = (float)',
                            source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                    }
                else
                    array_size.times {|i|
                        current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                            source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                    }
                end
            else
                type,array_name,array_size = return_type_name_size(subtext)

                if array_size == 1 && /double/ =~ type
                    current_output_file.print '    ',structure_name,'.',array_name,' = (double)',
                        source_structure_name,'[index]->',array_name,';'"\n"
                elsif array_size == 1
                    current_output_file.print '    ',structure_name,'.',array_name,' = ',
                        source_structure_name,'[index]->',array_name,';'"\n"
                elsif /double/ =~ type
                    array_size.times {|i|
                        current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = (double)',
                            source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                    }
                else
                    array_size.times {|i|
                        current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                            source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                    }
                end
            end

        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_copy_lean(structure_name,source_structure_name,current_read_file,current_output_file)
    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<LEAN>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)

            if array_size == 1 && /double/ =~ type
                current_output_file.print '    ',structure_name,'.',array_name,' = (float)',
                    source_structure_name,'[index]->',array_name,';'"\n"
            elsif array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    source_structure_name,'[index]->',array_name,';'"\n"
            elsif /double/ =~ type
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = (float)',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        source_structure_name,'[index]->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_copy_compact_elements(structure_name,pointer,source_structure_name,current_read_file,current_output_file)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    pointer,'->',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        pointer,'->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^\};/ =~ subtext
            break
        end
    end
end


def write_contains_of_structure_for_copy_compact_double_elements(structure_name,pointer,source_structure_name,current_read_file,current_output_file)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    pointer,'->',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        pointer,'->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_copy_compact_mix_elements(structure_name,pointer,source_structure_name,current_read_file,current_output_file)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<TMP>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    pointer,'->',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        pointer,'->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_contains_of_structure_for_copy_lean_elements(structure_name,pointer,source_structure_name,current_read_file,current_output_file)

    while subtext = current_read_file.gets do
        subtext.chomp!
        subtext.strip!
        subtext.gsub!('"','\"')

        if /^\/\// =~ subtext
            current_output_file.print "\n"
        elsif /(^\s*$|^\/\*)/ =~ subtext
            if /^\/\*/ =~ subtext
            else
                current_output_file.print "\n"
            end
        elsif /^#/ =~ subtext
            current_output_file.print subtext,"\n"
        elsif /<LEAN>/ =~ subtext
            current_output_file.print '// omit '+subtext+"\n"
        elsif /(bool|int|short|unsigned|float|double|void|struct|Struct)/ =~ subtext
            type,array_name,array_size = return_type_name_size(subtext)
            if array_size == 1
                current_output_file.print '    ',structure_name,'.',array_name,' = ',
                    pointer,'->',array_name,';'"\n"
            else
                array_size.times {|i|
                    current_output_file.print '    ',structure_name,'.',array_name,'[',i,']'+' = ',
                        pointer,'->',array_name,'[',i,']'';'"\n"
                }
            end
        elsif /^\};/ =~ subtext
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
        elsif /^\};/ =~ subtext
            break
        end
    end
end

def write_function_for_copy(function_name,source_structure_name,current_read_file,current_output_file,type)


        if type === "standard"
            current_output_file.print function_name,"\n"
            current_output_file.print "\n"

            split_structure_name = function_name.split(/\s+/)
            current_structure_name = split_structure_name[0]

            temporal_structuer_name = split_structure_name[0].sub("Struct","Temporarl")
            current_output_file.print "\t"+split_structure_name[0]+' '+temporal_structuer_name+";\n"
            current_output_file.print "\n"

            current_output_file.print "\t"+temporal_structuer_name+' = '+'*'+source_structure_name+'[index]'+"\n"
            current_output_file.print "\n"
            current_output_file.print "\treturn "+temporal_structuer_name+";\n"
            current_output_file.print "\n"
            current_output_file.print '}',"\n"
            current_output_file.print "\n"
        elsif type === "compact" || type === "compact_double" || type === "compact_mix" || type === "lean" 
            current_output_file.print function_name,"\n"
            current_output_file.print "\n"

            split_structure_name = function_name.split(/\s+/)
            current_structure_name = split_structure_name[1]

            temporal_structuer_name = split_structure_name[1].sub("Struct","Temp")
            current_output_file.print "\t"+split_structure_name[0]+' '+split_structure_name[1]+' '+temporal_structuer_name+";\n"
            current_output_file.print "\n"

            if type === "compact"
                write_contains_of_structure_for_copy_compact(temporal_structuer_name,
                    source_structure_name,current_read_file,current_output_file)
            elsif type === "compact_double"
                write_contains_of_structure_for_copy_compact_double(temporal_structuer_name,
                    source_structure_name,current_read_file,current_output_file)
            elsif type === "compact_mix"
                write_contains_of_structure_for_copy_compact_mix(temporal_structuer_name,
                    source_structure_name,current_read_file,current_output_file)
            elsif type === "lean"
                write_contains_of_structure_for_copy_lean(temporal_structuer_name,
                    source_structure_name,current_read_file,current_output_file)
            end
            current_output_file.print "\n"
            current_output_file.print "\treturn "+temporal_structuer_name+";\n"
            current_output_file.print "\n"
            current_output_file.print '}',"\n"
            current_output_file.print "\n"
        else
            STDOUT.print "File type error!!!\n"
        end

end

def write_function_for_copy_elements(function_name,variable,source_structure_name,current_read_file,current_output_file,type)

        current_output_file.print function_name+'Element'+'('+variable+'){'+"\n"
        current_output_file.print "\n"

        split_structure_name = function_name.split(/\s+/)
        current_structure_name = split_structure_name[1]

        temporal_structuer_name = split_structure_name[1].sub("Struct","Temp")
        current_output_file.print "\t"+split_structure_name[0]+' '+split_structure_name[1]+' '+temporal_structuer_name+";\n"
        current_output_file.print "\n"

        pointer_array = variable.split(/\s+/)
        pointer = pointer_array[2]

        if type === "compact"
            write_contains_of_structure_for_copy_compact_elements(temporal_structuer_name,pointer,
                    source_structure_name,current_read_file,current_output_file)
        elsif type === "compact_double"
            write_contains_of_structure_for_copy_compact_double_elements(temporal_structuer_name,pointer,
                    source_structure_name,current_read_file,current_output_file)
        elsif type === "compact_mix"
            write_contains_of_structure_for_copy_compact_mix_elements(temporal_structuer_name,pointer,
                    source_structure_name,current_read_file,current_output_file)
        elsif type === "lean"
            write_contains_of_structure_for_copy_lean_elements(temporal_structuer_name,pointer,
                    source_structure_name,current_read_file,current_output_file)
        end
        current_output_file.print "\n"
        current_output_file.print "\treturn "+temporal_structuer_name+";\n"
        current_output_file.print "\n"
        current_output_file.print '}',"\n"
        current_output_file.print "\n"
end

def read_function_from_copy(structure_name,current_read_file,current_output_file,current_header_file,version)

        type_of_structure = "Struct"+structure_name
        dest_structure_name = "Temporal"+structure_name
        source_structure_name = structure_name+"IO"
        
        if version === "standard"
            function_name = type_of_structure+" "+"CopyTemporalStructureTo"+
                structure_name+"("+type_of_structure+"IO"+" "+source_structure_name+"){"

            current_header_file.print type_of_structure+" "+"CopyTemporalStructureTo"+
                    structure_name+"("+type_of_structure+"IO"+" "+source_structure_name+");\n"

            current_output_file.print function_name,"\n"
            current_output_file.print "\n"

            current_output_file.print "\t"+type_of_structure+' '+dest_structure_name+";\n"
            current_output_file.print "\n"
        elsif version === "compact"
            function_name = type_of_structure+" "+"CopyTemporalStructureCompactTo"+
                structure_name+"Compact"+"(struct "+type_of_structure+"IO"+"Compact"+" "+source_structure_name+"Compact){"
            current_header_file.print type_of_structure+" "+"CopyTemporalStructureCompactTo"+
                structure_name+"Compact"+"(struct "+type_of_structure+"IO"+"Compact"+
                " "+source_structure_name+"Compact);"

            source_structure_name = source_structure_name+"Compact"

            current_output_file.print function_name,"\n"
            current_output_file.print "\n"

            current_output_file.print "\t"+type_of_structure+' '+dest_structure_name+";\n"
            current_output_file.print "\t"+'memset(&'+dest_structure_name+',0,sizeof('+type_of_structure+'))'+";\n"
            current_output_file.print "\n"

            read_contains_of_structure_from_copy(source_structure_name,
                    dest_structure_name,current_read_file,current_output_file,version)
        elsif version === "compact_double"
            function_name = type_of_structure+" "+"CopyTemporalStructureCompactDoubleTo"+
                structure_name+"CompactDouble"+"(struct "+type_of_structure+"IO"+"CompactDouble"+" "+source_structure_name+"CompactDouble){"
            current_header_file.print type_of_structure+" "+"CopyTemporalStructureCompactDoubleTo"+
                structure_name+"CompactDouble"+"(struct "+type_of_structure+"IO"+"CompactDouble"+
                " "+source_structure_name+"CompactDouble);"

            source_structure_name = source_structure_name+"CompactDouble"

            current_output_file.print function_name,"\n"
            current_output_file.print "\n"

            current_output_file.print "\t"+type_of_structure+' '+dest_structure_name+";\n"
            current_output_file.print "\t"+'memset(&'+dest_structure_name+',0,sizeof('+type_of_structure+'))'+";\n"
            current_output_file.print "\n"

            read_contains_of_structure_from_copy(source_structure_name,
                    dest_structure_name,current_read_file,current_output_file,version)
        elsif version === "compact_mix"
            function_name = type_of_structure+" "+"CopyTemporalStructureCompactMixTo"+
                structure_name+"CompactMix"+"(struct "+type_of_structure+"IO"+"CompactMix"+" "+source_structure_name+"CompactMix){"
            current_header_file.print type_of_structure+" "+"CopyTemporalStructureCompactMixTo"+
                structure_name+"CompactMix"+"(struct "+type_of_structure+"IO"+"CompactMix"+
                " "+source_structure_name+"CompactMix);"

            source_structure_name = source_structure_name+"CompactMix"

            current_output_file.print function_name,"\n"
            current_output_file.print "\n"

            current_output_file.print "\t"+type_of_structure+' '+dest_structure_name+";\n"
            current_output_file.print "\t"+'memset(&'+dest_structure_name+',0,sizeof('+type_of_structure+'))'+";\n"
            current_output_file.print "\n"

            read_contains_of_structure_from_copy(source_structure_name,
                    dest_structure_name,current_read_file,current_output_file,version)
        elsif version === "lean"
            function_name = type_of_structure+" "+"CopyTemporalStructureLeanTo"+
                structure_name+"Lean"+"(struct "+type_of_structure+"IO"+"Lean"+" "+source_structure_name+"Lean){"
            current_header_file.print type_of_structure+" "+"CopyTemporalStructureLeanTo"+
                structure_name+"Lean"+"(struct "+type_of_structure+"IO"+"Lean"+
                " "+source_structure_name+"Lean);"

            source_structure_name = source_structure_name+"Lean"

            current_output_file.print function_name,"\n"
            current_output_file.print "\n"

            current_output_file.print "\t"+type_of_structure+' '+dest_structure_name+";\n"
            current_output_file.print "\t"+'memset(&'+dest_structure_name+',0,sizeof('+type_of_structure+'))'+";\n"
            current_output_file.print "\n"

            read_contains_of_structure_from_copy_lean(source_structure_name,
                    dest_structure_name,current_read_file,current_output_file,version)
        else
            STDOUT.print "File type error!!\n"
        end

        current_output_file.print "\n"
        current_output_file.print "\treturn "+dest_structure_name+";\n"

        current_output_file.print "\n"
        current_output_file.print '}',"\n"
        current_output_file.print "\n"

end

=begin
First, this script writes "Standard Output" members in structures.
Next, the script writes "Compact Output" members in the structures.
At last, generate functions, which is copying from ordinary structure to temporal structures.
=end

$stderr.print "Make StructuresForIO.h/c, Start!\n"
read_file = open(filename)
output_header = open("StructuresForIO.h","w")
output_header.print "\n"
output_header.print '// This file is generated at ',Time.now,"\n"
output_header.print "\n"

output_body = open("StructuresForIO.c","w")
output_body.print "#include \"config.h\"\n"
output_body.print "#include \"DataStructures.h\"\n"
output_body.print "#include \"StructuresForIO.h\"\n"
output_body.print "\n"
output_body.print '// This file was generated at ',Time.now,"\n"
output_body.print "\n"

while text = read_file.gets do
    text.chomp!
    text.gsub!('"','\"')
    case text
    when "struct StructPbody_tag{"

        filepointer = read_file.pos

        # Compact format
        write_definition_of_structure("struct StructPbodyIOCompact{",read_file,output_header,"compact")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPbodyIOCompact CopyPbodyToTemporalStructureCompact(const int index){",
                "Pbody",read_file,output_body,"compact")
        output_header.print "struct StructPbodyIOCompact CopyPbodyToTemporalStructureCompact(const int index);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pbody",read_file,output_body,output_header,"compact")
        output_header.print "\n\n"


        # Compact_double format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPbodyIOCompactDouble{",read_file,output_header,"compact_double")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPbodyIOCompactDouble CopyPbodyToTemporalStructureCompactDouble(const int index){",
                "Pbody",read_file,output_body,"compact_double")
        output_header.print "struct StructPbodyIOCompactDouble CopyPbodyToTemporalStructureCompactDouble(const int index);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pbody",read_file,output_body,output_header,"compact_double")
        output_header.print "\n\n"

        # Compact_mix format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPbodyIOCompactMix{",read_file,output_header,"compact_mix")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPbodyIOCompactMix CopyPbodyToTemporalStructureCompactMix(const int index){",
                "Pbody",read_file,output_body,"compact_mix")
        output_header.print "struct StructPbodyIOCompactMix CopyPbodyToTemporalStructureCompactMix(const int index);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pbody",read_file,output_body,output_header,"compact_mix")
        output_header.print "\n\n"

        # Lean format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPbodyIOLean{",read_file,output_header,"lean")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPbodyIOLean CopyPbodyToTemporalStructureLean(const int index){",
                "Pbody",read_file,output_body,"lean")
        output_header.print "struct StructPbodyIOLean CopyPbodyToTemporalStructureLean(const int index);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pbody",read_file,output_body,output_header,"lean")
        output_header.print "\n\n"

    when "struct StructPhydro_tag{"

        filepointer = read_file.pos

        # Compact format
        write_definition_of_structure("struct StructPhydroIOCompact{",read_file,output_header,"compact")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompact(const int index){",
                "Phydro",read_file,output_body,"compact")
        output_header.print "struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompact(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompact",
                "StructPhydroptr const Ph","Pbody",read_file,output_body,"compact")
        output_header.print "struct StructPhydroIOCompact CopyPhydroToTemporalStructureCompactElement(StructPhydroptr const Ph);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Phydro",read_file,output_body,output_header,"compact")
        output_header.print "\n\n"


        # Compact_double format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPhydroIOCompactDouble{",read_file,output_header,"compact_double")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPhydroIOCompactDouble CopyPhydroToTemporalStructureCompactDouble(const int index){",
                "Phydro",read_file,output_body,"compact_double")
        output_header.print "struct StructPhydroIOCompactDouble CopyPhydroToTemporalStructureCompactDouble(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPhydroIOCompactDouble CopyPhydroToTemporalStructureCompactDouble",
                "StructPhydroptr const Ph","Pbody",read_file,output_body,"compact_double")
        output_header.print "struct StructPhydroIOCompactDouble CopyPhydroToTemporalStructureCompactDoubleElement(StructPhydroptr const Ph);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Phydro",read_file,output_body,output_header,"compact_double")
        output_header.print "\n\n"


        # Compact_mix format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPhydroIOCompactMix{",read_file,output_header,"compact_mix")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPhydroIOCompactMix CopyPhydroToTemporalStructureCompactMix(const int index){",
                "Phydro",read_file,output_body,"compact_mix")
        output_header.print "struct StructPhydroIOCompactMix CopyPhydroToTemporalStructureCompactMix(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPhydroIOCompactMix CopyPhydroToTemporalStructureCompactMix",
                "StructPhydroptr const Ph","Pbody",read_file,output_body,"compact_mix")
        output_header.print "struct StructPhydroIOCompactMix CopyPhydroToTemporalStructureCompactMixElement(StructPhydroptr const Ph);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Phydro",read_file,output_body,output_header,"compact_mix")
        output_header.print "\n\n"


        # Lean format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPhydroIOLean{",read_file,output_header,"lean")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPhydroIOLean CopyPhydroToTemporalStructureLean(const int index){",
                "Phydro",read_file,output_body,"lean")
        output_header.print "struct StructPhydroIOLean CopyPhydroToTemporalStructureLean(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPhydroIOLean CopyPhydroToTemporalStructureLean",
                "StructPhydroptr const Ph","Pbody",read_file,output_body,"lean")
        output_header.print "struct StructPhydroIOLean CopyPhydroToTemporalStructureLeanElement(StructPhydroptr const Ph);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Phydro",read_file,output_body,output_header,"lean")
        output_header.print "\n\n"

    when "struct StructPstar_tag{"

        filepointer = read_file.pos

        # Compact format
        write_definition_of_structure("struct StructPstarIOCompact{",read_file,output_header,"compact")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPstarIOCompact CopyPstarToTemporalStructureCompact(const int index){",
                "Pstar",read_file,output_body,"compact")
        output_header.print "struct StructPstarIOCompact CopyPstarToTemporalStructureCompact(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPstarIOCompact CopyPstarToTemporalStructureCompact",
                "StructPstarptr const Ps","Pbody",read_file,output_body,"compact")
        output_header.print "struct StructPstarIOCompact CopyPstarToTemporalStructureCompactElement(StructPstarptr const Ps);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pstar",read_file,output_body,output_header,"compact")
        output_header.print "\n\n"


        # Compact_double format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPstarIOCompactDouble{",read_file,output_header,"compact_double")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPstarIOCompactDouble CopyPstarToTemporalStructureCompactDouble(const int index){",
                "Pstar",read_file,output_body,"compact_double")
        output_header.print "struct StructPstarIOCompactDouble CopyPstarToTemporalStructureCompactDouble(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPstarIOCompactDouble CopyPstarToTemporalStructureCompactDouble",
                "StructPstarptr const Ps","Pbody",read_file,output_body,"compact_double")
        output_header.print "struct StructPstarIOCompactDouble CopyPstarToTemporalStructureCompactDoubleElement(StructPstarptr const Ps);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pstar",read_file,output_body,output_header,"compact_double")
        output_header.print "\n\n"


        # Compact_mix format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPstarIOCompactMix{",read_file,output_header,"compact_mix")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPstarIOCompactMix CopyPstarToTemporalStructureCompactMix(const int index){",
                "Pstar",read_file,output_body,"compact_mix")
        output_header.print "struct StructPstarIOCompactMix CopyPstarToTemporalStructureCompactMix(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPstarIOCompactMix CopyPstarToTemporalStructureCompactMix",
                "StructPstarptr const Ps","Pbody",read_file,output_body,"compact_mix")
        output_header.print "struct StructPstarIOCompactMix CopyPstarToTemporalStructureCompactMixElement(StructPstarptr const Ps);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pstar",read_file,output_body,output_header,"compact_mix")
        output_header.print "\n\n"


        # Lean format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPstarIOLean{",read_file,output_header,"lean")

        #cp function

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPstarIOLean CopyPstarToTemporalStructureLean(const int index){",
                "Pstar",read_file,output_body,"lean")
        output_header.print "struct StructPstarIOLean CopyPstarToTemporalStructureLean(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPstarIOLean CopyPstarToTemporalStructureLean",
                "StructPstarptr const Ps","Pbody",read_file,output_body,"lean")
        output_header.print "struct StructPstarIOLean CopyPstarToTemporalStructureLeanElement(StructPstarptr const Ps);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Pstar",read_file,output_body,output_header,"lean")
        output_header.print "\n\n"

    when "struct StructPsink_tag{"

        filepointer = read_file.pos

        # Compact format
        write_definition_of_structure("struct StructPsinkIOCompact{",read_file,output_header,"compact")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompact(const int index){",
                "Psink",read_file,output_body,"compact")
        output_header.print "struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompact(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompact",
                "StructPsinkptr const Psk","Pbody",read_file,output_body,"compact")
        output_header.print "struct StructPsinkIOCompact CopyPsinkToTemporalStructureCompactElement(StructPsinkptr const Psk);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Psink",read_file,output_body,output_header,"compact")
        output_header.print "\n\n"


        # Compact_double format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPsinkIOCompactDouble{",read_file,output_header,"compact_double")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPsinkIOCompactDouble CopyPsinkToTemporalStructureCompactDouble(const int index){",
                "Psink",read_file,output_body,"compact_double")
        output_header.print "struct StructPsinkIOCompactDouble CopyPsinkToTemporalStructureCompactDouble(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPsinkIOCompactDouble CopyPsinkToTemporalStructureCompactDouble",
                "StructPsinkptr const Psk","Pbody",read_file,output_body,"compact_double")
        output_header.print "struct StructPsinkIOCompactDouble CopyPsinkToTemporalStructureCompactDoubleElement(StructPsinkptr const Psk);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Psink",read_file,output_body,output_header,"compact_double")
        output_header.print "\n\n"


        # Compact_mix format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPsinkIOCompactMix{",read_file,output_header,"compact_mix")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPsinkIOCompactMix CopyPsinkToTemporalStructureCompactMix(const int index){",
                "Psink",read_file,output_body,"compact_mix")
        output_header.print "struct StructPsinkIOCompactMix CopyPsinkToTemporalStructureCompactMix(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPsinkIOCompactMix CopyPsinkToTemporalStructureCompactMix",
                "StructPsinkptr const Psk","Pbody",read_file,output_body,"compact_mix")
        output_header.print "struct StructPsinkIOCompactMix CopyPsinkToTemporalStructureCompactMixElement(StructPsinkptr const Psk);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Psink",read_file,output_body,output_header,"compact_mix")
        output_header.print "\n\n"


        # Lean format
        read_file.seek(filepointer,IO::SEEK_SET)
        write_definition_of_structure("struct StructPsinkIOLean{",read_file,output_header,"lean")

        #cp function
        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy("struct StructPsinkIOLean CopyPsinkToTemporalStructureLean(const int index){",
                "Psink",read_file,output_body,"lean")
        output_header.print "struct StructPsinkIOLean CopyPsinkToTemporalStructureLean(const int index);\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        write_function_for_copy_elements("struct StructPsinkIOLean CopyPsinkToTemporalStructureLean",
                "StructPsinkptr const Psk","Pbody",read_file,output_body,"lean")
        output_header.print "struct StructPsinkIOLean CopyPsinkToTemporalStructureLeanElement(StructPsinkptr const Psk);\n"
        output_header.print "\n"

        read_file.seek(filepointer,IO::SEEK_SET)
        read_function_from_copy("Psink",read_file,output_body,output_header,"lean")
        output_header.print "\n\n"

    end
end

output_header.print "\n"
output_body.print "\n"

read_file.close 
output_header.close 
output_body.close 

$stderr.print "Make StructIOAG.h, Finished!\n"


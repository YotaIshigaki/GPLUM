def getName(val):
    if(val.split()[0] == "static"):
        return val.split()[2]
    else:
        return val.split()[1]
    
def getType(val):
    if(val.split()[0] == "static"):
        return val.split()[1]
    else:
        return val.split()[0]
    
def getStrCopyFromFp(fp, interaction, is_epi):
    keyword = "EPIVALS"
    if(not is_epi):
        keyword = "EPJVALS"
    str_copy_from_fp = ""
    for i in fp:
        for j in range(len(i["INTERACTION"])):
            if( i["INTERACTION"][j] == interaction["NAME"]):
                str_copy_from_fp += "    void copyFromFP (const FP" + i["NAME"] + " & x) {\n"
                for k in interaction[keyword]:
                    if(k.split()[0] == "static"):
                        continue
                    val_name = getName(k)
                    str_copy_from_fp += "        this->" + val_name + " = x." + val_name + ";\n"
                str_copy_from_fp += "    }\n"
    return str_copy_from_fp


def getStrAccessor(str_val):
    str_accessor = ""
    if(str_val == ""):
        return str_accessor
    keyword = "Pos"
    if( getName(str_val) == "r_search"):
        keyword = "RSearch"
    elif( getName(str_val) == "mass"):
        keyword = "Charge"
    if(getName(str_val) != ""):
        type_val = getType(str_val)
        name_val = getName(str_val)
        str_accessor = """
    {TYPEVAL} get{KEYWORD} () const {{
        return {NAMEVAL};
    }}
        """.format(TYPEVAL=type_val, NAMEVAL=name_val, KEYWORD=keyword)
        if(getName(str_val) == "pos"):
                    str_accessor += """
    void set{KEYWORD} (const {TYPEVAL} & x) {{
        {NAMEVAL} = x;
    }}
        """.format(TYPEVAL=type_val, NAMEVAL=name_val, KEYWORD=keyword)
    return str_accessor
            
class Force:
    def __init__(self, interaction):
        self.interaction = interaction
    def dump(self):
        str_force_vals = ""
        for i in self.interaction["FORCEVALS"]:
            str_force_vals += "    " + i + ";\n"

        str_clear = "    void clear () {\n"
        for i in self.interaction["FORCEVALS"]:
            str_clear += "        this->" + i.split()[1] + " = 0.0; \n"
        str_clear += "    }\n"
            
        string = """
class FORCE{NAME} {{
public:
{FORCEVALS}
{CLEAR}
}};
        """.format(NAME = self.interaction["NAME"], FORCEVALS=str_force_vals, CLEAR=str_clear)
        print(string)
        
class Fp:
    def __init__(self, fp, interaction):
        self.fp = fp
        self.interaction = interaction
        self.epi_vals = []
        self.epj_vals = []
        self.force_vals = []
        for i in self.fp["INTERACTION"]:
            for j in self.interaction:
                if(i == j["NAME"]):
                    for k in j["EPIVALS"]:
                        self.epi_vals.append(k)
                    for k in j["EPJVALS"]:
                        self.epj_vals.append(k)
                    for k in j["FORCEVALS"]:
                        self.force_vals.append(k)
        self.ep_vals = list(set(self.epi_vals + self.epj_vals))
        
    def dump(self):
        str_pos     = ""
        str_mass     = ""
        str_r_search = ""
        str_othervals = ""
        for i in self.fp["OTHERVALS"]:
            str_othervals += "    " + i + ";\n"
            if(getName(i) == "pos"):
                str_pos = i
            if(getName(i) == "mass"):
                str_mass = i
            if(getName(i) == "r_search"):
                str_r_search = i
            
        str_ep_vals = ""
        for i in self.ep_vals:
            str_ep_vals += "    " + i + ";\n"
            if(getName(i) == "pos"):
                str_pos = i
            if(getName(i) == "mass"):
                str_mass = i
            if(getName(i) == "r_search"):
                str_r_search = i
                
        str_force_vals = ""
        for i in self.force_vals:
            str_force_vals += "    " + i + ";\n"
            if(getName(i) == "pos"):
                str_pos = i
            if(getName(i) == "mass"):
                str_mass = i
            if(getName(i) == "r_search"):
                str_r_search = i

        str_accessor = getStrAccessor(str_pos)
        str_accessor += getStrAccessor(str_mass)
        str_accessor += getStrAccessor(str_r_search)

        str_copy_from_force = ""
        for i in self.interaction:
            for j in self.fp["INTERACTION"]:
                if( i["NAME"] == j):
                    str_copy_from_force += "    void copyFromForce (const FORCE" + i["NAME"] + " & x) {\n"
                    #for k in self.force_vals:
                    for k in i["FORCEVALS"]:
                        str_copy_from_force += "        this->" + getName(k) + " = x." + getName(k) + ";\n"
                    str_copy_from_force += "    }\n"
            
        string = """
class FP{NAME} {{
public:
{EPVALS}
{FORCEVALS}
{OTHERVALS}
{ACCESSOR}
{COPYFROMFORCE}
}};
        """.format(NAME = self.fp["NAME"], EPVALS = str_ep_vals, FORCEVALS=str_force_vals, OTHERVALS = str_othervals, ACCESSOR=str_accessor, COPYFROMFORCE = str_copy_from_force)
        print(string)


class Ep:
    def __init__(self, interaction, fp, is_epi):
        self.interaction = interaction
        self.fp = fp
        self.is_epi = is_epi
    def dump(self):
        str_vals = "EPIVALS"
        name_type = "EPI" + self.interaction["NAME"]
        if(not self.is_epi):
            str_vals = "EPJVALS"
            name_type = "EPJ" + self.interaction["NAME"]
            
        str_epi_vals = ""
        str_mass = ""
        str_pos = ""
        str_r_search = ""
        for i in self.interaction[str_vals]:
            str_epi_vals += "    " + i + ";\n"
            if(getName(i) == "pos"):
                str_pos = i
            if(getName(i) == "mass"):
                str_mass = i
            if(getName(i) == "r_search"):
                str_r_search = i
                
        str_accessor = getStrAccessor(str_pos)
        str_accessor += getStrAccessor(str_mass)
        str_accessor += getStrAccessor(str_r_search)
            
        is_epi = True            
        str_copy_from_fp = getStrCopyFromFp(self.fp, self.interaction, is_epi)
        
        string = """
class {NAMETYPE} {{
public:
{EPIVALS}
{ACCESSOR}
{COPYFROMFP}
}};
        """.format(NAMETYPE = name_type, EPIVALS=str_epi_vals, ACCESSOR=str_accessor, COPYFROMFP=str_copy_from_fp)
        print(string)
        
class ForceFunc:
    def __init__(self, interaction):
        self.interaction = interaction
    def dump(self):
        name_interaction = self.interaction["NAME"]
        name_epi   = "EPI"   + name_interaction
        name_epj   = "EPJ"   + name_interaction
        name_force = "FORCE" + name_interaction
        string = """
void Calc{INTERACTION}EpEp ( 
    const {EPINAME} * epi,
    const PS::S32 n_ip,
    const {EPJNAME} * epj,
    const PS::S32 n_jp,
    {FORCENAME} * force){{
    for(PS::S32 i = 0; i < n_ip; i++){{
        for(PS::S32 j = 0; j < n_jp; j++){{
        // do something
        }}
    }}
}}

template<class Tspj>
void Calc{INTERACTION}EpSp ( 
    const {EPINAME} * epi,
    const PS::S32 n_ip,
    const Tspj * spj,
    const PS::S32 n_jp,
    {FORCENAME} * force){{
    for(PS::S32 i = 0; i < n_ip; i++){{
        for(PS::S32 j = 0; j < n_jp; j++){{
        // do something
        }}
    }}
}}
        """.format(INTERACTION=name_interaction, EPINAME=name_epi, EPJNAME=name_epj, FORCENAME=name_force)
        print(string)

        
import sys
import yaml

file_name = sys.argv[1]
print(file_name)
file = open(file_name)

data = yaml.load(file)

#print(data)
n_interaction = len(data["INTERACTION"])
n_fp = len(data["FP"])

force = []
for i in range(n_interaction):
    force.append( Force(data["INTERACTION"][i]) )
    force[i].dump()

fp = []
for i in range(n_fp):
    fp.append( Fp(data["FP"][i], data["INTERACTION"]) )
    fp[i].dump()

epi = []    
for i in range(n_interaction):
    is_epi = True
    epi.append( Ep(data["INTERACTION"][i], data["FP"], is_epi) )
    epi[i].dump()

epj = []    
for i in range(n_interaction):
    is_epi = False
    epj.append( Ep(data["INTERACTION"][i], data["FP"], is_epi) )
    epj[i].dump()


force_func = []
for i in range(n_interaction):
    force_func.append( ForceFunc(data["INTERACTION"][i]) )
    force_func[i].dump()
    


import re

Obsidian = open("data.json","rt")
Latex = open("latex.json","wt")
Obmod = Obsidian.readlines()[18]
Obmod = Obmod.replace("\"customShorthand_parameter\":", "")
Obmod = Obmod.replace("\"","")
Obmod = Obmod.replace(" ","")
Obmod = Obmod.replace("#cursor","$1")
Obmod = Obmod.replace("#tab","$2")
Obmod = Obmod.replace("\n","")
Obsplit = Obmod.split(";")
Obsplit = Obsplit[:-1]
Latex.write("{\n")
for x in Obsplit:
    x = x.replace("\\n","")
    testiculus = x.split(":")
    
    result = "\t\"" + testiculus[0] +"\":{\n\t\t\"prefix\": \""+testiculus[0] +"\",\n\t\t\"body\": [\""+testiculus[1]+ "\"],\n\t\t\"description\": \"\"\n\t}\n"
    Latex.write(result)
Latex.write("}")
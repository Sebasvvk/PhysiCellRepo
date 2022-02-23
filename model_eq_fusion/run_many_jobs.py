# This script provides simple parameter exploration functionality. The script creates
# a new folder (subdirectory) for each set of parameters, makes changes to a default 
# configuration (.xml) file using specified parameter values (in an accompanying .txt file),
# copies the new config file into the new folder, then
# runs the simulation (in the background) which writes results into the new folder.
# 

import xml.etree.ElementTree as ET
from shutil import copyfile
import os
import sys
import subprocess



exec_pgm='./project'
background_str = " &"  # works on Unix
background_str=" "
if sys.platform == 'win32':
    background_str = ""


xml_file_in = 'config/PhysiCell_settings.xml'
xml_file_tmp = 'config/tmp.xml'
copyfile(xml_file_in, xml_file_tmp)

first_time = True
output_dirs = []


adhesions=[0.4,1.2,4]
speeds=[0,0.1,0.2,0.5]
#follower_speeds=[0.1,0.2,0.3,0.4]
for adhesion in adhesions:
    for speed in speeds:
        pre='onesphere_adh_'+str(adhesion)+'_speed_'+str(speed)
        ### make folder ###
        folder_name='output_'+pre
        output_dirs.append(folder_name)
        if (not os.path.exists(folder_name)):
            print("--- parsed 'folder', makedir " + folder_name)
            os.makedirs(folder_name)
        ### end make folder

        
        ### replace parameters ###
        tree = ET.parse(xml_file_tmp)
        xml_root = tree.getroot()
        xml_root.find('.//' + 'folder').text = folder_name ### don't forget this to save data to the correct folder
        #xml_root.find('.//' + 'r_al0').text = str(r_al0)
        xml_root.find('.//' + 'cell_definition[@name="default"]').find('./phenotype/mechanics/cell_cell_adhesion_strength').text=str(adhesion)
        xml_root.find('.//' + 'cell_definition[@name="second"]').find('./phenotype/mechanics/cell_cell_adhesion_strength').text=str(adhesion)
        xml_root.find('.//' + 'cell_definition[@name="default"]').find('./phenotype/motility/speed').text=str(speed)
        xml_root.find('.//' + 'cell_definition[@name="second"]').find('./phenotype/motility/speed').text=str(speed)

        #xml_root.find('.//' + 'cell_definition[@name="follower cell"]').find('./phenotype/motility/speed').text=str(follower_speed)
        

        ### end replace parameters###

        ### write config file ###
        xml_file_out = os.path.join(folder_name, 'config.xml')  # copy config file into the output dir
        print('---write config file (and start sim): ', xml_file_out)
        tree.write(xml_file_out)   # will create folder_name/config.xml
        log_file = folder_name + ".log"  
        cmd =  exec_pgm + " " + xml_file_out + " > " + log_file + " " + background_str
        print("----- cmd = ",cmd)
        # os.system(cmd)   # <------ Execute the simulation
        # subprocess.Popen([exec_pgm, xml_file_out])
        with open(log_file,"w") as outf:
            subprocess.run([exec_pgm, xml_file_out],stdout=outf)
        cmd="sh makemovie.sh "+folder_name
#            log_file = folder_name + ".logmovie"  
#            with open(log_file,"w") as outf:
#                subprocess.run(['sh', 'makemovie.sh '+folder_name],stdout=outf)
        os.system(cmd)
        cmd="python3 to_dump.py "+folder_name
        os.system(cmd)



print("\n ------\n Your output results will appear in these directories:\n   ",output_dirs)
print("and check for a .log file of each name for your terminal output from each simulation.\n")

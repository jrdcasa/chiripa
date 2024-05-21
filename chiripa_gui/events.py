import PySimpleGUI as sg
import webbrowser
import datetime
import os
import ast
import re
from menu_layout import menu_layout
import clipboard
import subprocess
import glob
import time
import chiripa as chi

# DICTIONARY
dict_topology = {'-TOPOCHECKBOX1-': ['-TOPOFILE1-', '-BROWSETOPO1-' ],
                 '-TOPOCHECKBOX2-': ['-TOPOFILE2-', '-BROWSETOPO2-' ],}
dict_coord = {'-COORDFILE1-': ['1', '-NAME1-', '-TOPOCHECKBOX1-', '-TOPOFILE1-'],
              '-COORDFILE2-': ['2', '-NAME2-', '-TOPOCHECKBOX2-', '-TOPOFILE2-']}
d2 = {'-TOPOFILE1-': ['1'], '-TOPOFILE2-': ['2']}
dict_button_script = {'-BUTTONCREATESCRIPT-': ['-NAME1-', '-NAME2-',
                                               '-COORDFILE1-', '-COORDFILE2-',
                                               '-TOPOFILE1-', '-TOPOFILE2-',
                                               '-QMPACK-', '-QMPATHEXE-', '-QMPATHSCRATCH-',
                                               '-QMMETHOD-','-QMBASISSET-', '-QMTASK-', '-SERVERNAME-',
                                               '-SERVERREMOTEDIR-', '-SERVERLOCALDIR-'],
                      '-BUTTONRUN-'         :  ['-NAME1-', '-NAME2-',
                                               '-COORDFILE1-', '-COORDFILE2-',
                                               '-TOPOFILE1-', '-TOPOFILE2-',
                                               '-QMPACK-', '-QMPATHEXE-', '-QMPATHSCRATCH-',
                                               '-QMMETHOD-','-QMBASISSET-', '-QMTASK-', '-SERVERNAME-',
                                               '-SERVERREMOTEDIR-', '-SERVERLOCALDIR-'],
                      }

dict_Zsamples = {'-ZCheckBox-' : ['-ZSAMPLES-', '-ZTRIALS-', '-ZDEBUG-', '-ZNONBONDED-']}
dict_QMoptions = {'dict_QMoptions-' : ['-QMPACK-', '-QMPATHEXE-', '-QMPATHSCRATCH-', '-QMCONF-',
                                       '-QMCHARGE-', '-QMMULTIPLICITY-', '-QMTASK-', '-QMMETHOD-',
                                       '-QMBASISSET-', '-QMOTHER-', '-QMMEM-'
                                       ]}
dict_Serveroptions = {'localhost' : ['-SERVERNAME-', '-SERVERREMOTEDIR-', '-SERVERLOCALDIR-'],
                      'slurm' : ['-SERVERNAME-', '-SERVERQUEUE-', '-SERVERUSER-', '-SERVERKEYFILE-',
                                 '-SERVERREMOTEDIR-', '-SERVERLOCALDIR-', '-SERVERPARTITION-',
                                 '-SERVERNODELIST-', '-SERVERPARTMASTER-', '-SERVERPARTMASTERNODE-']}

# =============================================================================
def molecule_web(window, molpath):

    """
    Create the index.html used by `ngl <http://nglviewer.org/#ngl>`_ viewer embebed  in a web to
    represent the structures. The index_template.html file must exist in the www directory.
    """

    # Open index_template.html

    string_files = []
    ext_files = []
    try:
        for ifile in range(2):
            with open(molpath[ifile], 'r') as f:
                l = f.read()
                string_files.append(repr(l))
                ext_files.append(molpath[ifile].split(".")[-1])

        with open("./www/index_template.html", 'r') as f:

            lines = f.read()
            lines = lines.replace("#FILE1STRING#", string_files[0])
            lines = lines.replace("#FILE2STRING#", string_files[1])
            lines = lines.replace("#FILE1EXT#", ext_files[0])
            lines = lines.replace("#FILE2EXT#", ext_files[1])

        # Save file
        with open("./www/index.html", 'w') as f:
            f.writelines(lines)

        webbrowser.open("./www/index.html", new=1)

    except FileNotFoundError:

        loc_main = window.current_location()
        size = window.size
        loc_center = (loc_main[0] + size[0] / 2, loc_main[1] + size[1] / 2)
        sg.popup("File {} does not exist".format(molpath[ifile]), title="Warning", location = loc_center)

# =============================================================================
def enable_disable_topology(window, dict_list):

    """ Enable/Disable Inputs and buttons as function of topology checkbox"""

    for key, elem_list in dict_list.items():
        if window[key].get():
            isdisabled = True
        else:
            isdisabled = False
        for item in elem_list:
            window[item].update(disabled=isdisabled)

# =============================================================================
def enable_disable_buttonscript(window, dict_list, dict_Zsamples, dict_Serveroptions):

    """ Enable/Disable Python Script Button and run button"""
    isdisabled = True
    l = []
    for key, elem_list in dict_list.items():
        for item in elem_list:
            if window[item].get() == '':
                l.append("True")
            else:
                l.append("False")

        for key2, elem_list2 in dict_Zsamples.items():
            if window['-ZCheckBox-'].get():
                for item in elem_list2:
                    if window[item].get() == '':
                        l.append("True")
                    else:
                        l.append("False")

        for key3, elem_list3 in dict_QMoptions.items():
            if window['-IECheckBox-'].get():
                for item in elem_list3:
                    if item == '-QMOTHER-' or item == '-QMMEM-' :
                        continue
                    if window[item].get() == '':
                        l.append("True")
                    else:
                        l.append("False")

        for key4, elem_list4 in dict_Serveroptions.items():
            if key4 == "localhost" and window["-RADIO1a-"].get():
                for item in elem_list4:
                    if window[item].get() == '':
                        l.append("True")
                    else:
                        l.append("False")
            elif key4 == "slurm" and window["-RADIO1b-"].get():
                for item in elem_list4:
                    if window[item].get() == '':
                        l.append("True")
                    else:
                        l.append("False")
            else:
                pass

        if all( i=="False" for i in l):
            isdisabled = False

        window[key].update(disabled=isdisabled)
        
        # Unable/Disable Menu/Export...
        if not isdisabled:
            menu_def = [['File', ['Import...', 'Export...', 'Exit']],
                        ['Edit', ['Clean Form'], ],
                        ['Help', 'About...'], ]
            window["-MENU-"].update(menu_definition=menu_def)
        else:
            menu_def = [['File', ['Import...', '!Export...', 'Exit']],
                        ['Edit', ['Clean Form'], ],
                        ['Help', 'About...'], ]
            window["-MENU-"].update(menu_definition=menu_def)

    return None

# =============================================================================
def enable_disable_Zoptions(window, dict_Zsamples):

    for key, values in dict_Zsamples.items():
            for item in values:
                if window["-ZCheckBox-"].get():
                    window[item].update(disabled=False)
                else:
                    window[item].update(disabled=True)

# =============================================================================
def create_script_chiripa(window, save=False, filename="script_chiripa.py"):

    l  = "# Import chiripa package\n"
    l += "import chiripa as chi\n"
    l += "\n"
    l += "# Keyword dictionary\n"
    l += "inputdict = {\n"
    l += "      \'names\'                  : [\'{}\', \'{}\'],\n".format(window['-NAME1-'].get(),window['-NAME2-'].get())
    l += "      \'filecoords\'             : [\'{}\', \'{}\'],\n".format(window['-COORDFILE1-'].get(),window['-COORDFILE2-'].get())
    l += "      \'filetop\'                : [\'{}\', \'{}\'],\n".format(window['-TOPOFILE1-'].get(),window['-TOPOFILE2-'].get())
    l += "      \'coordination_numbers_Z\' : {},\n".format("True" if window['-ZCheckBox-'].get() != 0 else "False")
    l += "      \'calculate_volume\'       : {},\n".format("True" if window['-VolCheckBox-'].get() != 0 else "False")
    l += "      \'interaction_energy\'     : {},\n".format("True" if window['-IECheckBox-'].get() != 0 else "False")
    l += "      \'montecarlo_run\'         : {},\n".format("True" if window['-MCCheckBox-'].get() != 0 else "False")
    l += "      \'boxl\'                   : {},\n".format(window['-BOXL-'].get())
    if window['-ZCheckBox-'].get() == 1:
        l += "      \'Z_parameters\'           : {{ \'Z_samples\'            : {},\n".format(window['-ZSAMPLES-'].get())
        l += "                                   \'Z_puttrialmonomers\'   : {},\n".format(window['-ZTRIALS-'].get())
        l += "                                   \'Z_debug\'              : {},\n".format("True" if window['-ZDEBUG-'].get() != 0 else "False")
        l += "                                   \'Z_nonbonded\'          : \"{}\" }},\n".format(window['-ZNONBONDED-'].get())


    if window['-MCCheckBox-'].get() == 1:
        l += "      \'montecarlo_parameters\'  :  {{ \'temperatures\'            : [{}, {}, {}],\n".\
            format(window['-MCINITEMP-'].get(),window['-MCENDTEMP-'].get(), window['-MCSTEPTEMP-'].get() )
        l += "                                   \'mode\'   : \"{}\",\n".format(window['-MCMODE-'].get())
        l += "                                   \'pairs\'              : \"{}\",}},\n".format(window['-MCPAIRS-'].get())



    if window['-IECheckBox-'].get() == 1:
        l += "      \'energy_parameters\'      : {{ \'qm_engine\'            : \"{}\",\n".format(window['-QMPACK-'].get())
        l += "                                   \'qm_path_exe\'          : \"{}\",\n".format(window['-QMPATHEXE-'].get())
        l += "                                   \'qm_scratch_dir\'       : \"{}\",\n".format(window['-QMPATHSCRATCH-'].get())
        l += "                                   \'qm_charge\'            : {},\n".format(window['-QMCHARGE-'].get())
        l += "                                   \'qm_multiplicity\'      : {},\n".format(window['-QMMULTIPLICITY-'].get())
        l += "                                   \'qm_basisset\'          : \"{}\",\n".format(window['-QMBASISSET-'].get())
        l += "                                   \'qm_method\'            : \"{}\",\n".format(window['-QMMETHOD-'].get())
        l += "                                   \'qm_task\'              : \"{}\",\n".format(window['-QMTASK-'].get())
        if len(window['-QMMEM-'].get()) != 0:
            l += "                                   \'qm_memory_mb\'         : {},\n".format(window['-QMMEM-'].get())
        l += "                                   \'number_configurations\': {} }},\n".format(window['-QMCONF-'].get())

    l += "      \'server\'                 : {{ \'name\'          : \"{}\",\n".format(window['-SERVERNAME-'].get())
    l += "                                      \'remote_dir\'  : \"{}\",\n".format(window['-SERVERREMOTEDIR-'].get())
    l += "                                      \'local_dir\'   : \"{}\",\n".format(window['-SERVERLOCALDIR-'].get())
    l += "                                      \'ncpus\'       : \"{}\",\n".format(window['-SERVERCPU-'].get())

    if window['-RADIO1b-'].get():
        l += "                                      \'queue_system\'       : \"{}\",\n".format(window['-SERVERQUEUE-'].get())
        l += "                                      \'username\'           : \"{}\",\n".format(window['-SERVERUSER-'].get())
        l += "                                      \'username\'           : \"{}\",\n".format(window['-SERVERUSER-'].get())
        l += "                                      \'partition\'          : \"{}\",\n".format(window['-SERVERPARTITION-'].get())
        l += "                                      \'nodelist\'           : {},\n".format(window['-SERVERNODELIST-'].get().split(","))
        l += "                                      \'partitionmaster\'    : \"{}\",\n".format(window['-SERVERPARTMASTER-'].get())
        l += "                                      \'nodelistmaster\'     : \"{}\",}}\n".format(window['-SERVERPARTMASTERNODE-'].get())

    l += "}}\n"

    l += "c = chi.Chi_Universe(inputdict, output_screen=True)\n"

    if save and not filename is None:
        try:
            with open(filename, 'w') as f:
                f.writelines(l)
        except FileNotFoundError:
            pass

    return l

# =============================================================================
def import_file_to_gui(window, filename):

    """
    Use a python script for Chiripa and load the parameters in the GUI

    Args:
        filename: File containing chiripa  dictionary

    Returns:

    """

    # If filename None, the CANCEL button has been pushed
    if filename is None or len(filename) == 0:
        return False

    # Open file containing inputdict
    with open(filename, 'r') as fin:
        # lines = fin.readlines()
        lines = fin.read()

    working_dir, file = os.path.split(filename)
    prev_dir = os.path.abspath(os.path.join(os.path.dirname(filename), '..', ''))

    data = re.findall(r'{(?:.|\s)+}', lines)
    inputdict = ast.literal_eval(data[0].strip())

    # Read keys and fill up the forms in the gui
    for key, item in inputdict.items():
        if key == 'names':
            if isinstance(item, list) and len(item) == 2:
                window['-NAME1-'].update(item[0])
                window['-NAME2-'].update(item[1])
        elif key == 'filecoords':
            if isinstance(item, list) and len(item) == 2:
                if item[0].find(".") != 1:
                    f1 = item[0].replace(".",working_dir)
                if item[1].find(".") != 1:
                    f2 = item[1].replace(".",working_dir)
                if item[0].find("..") != 1:
                    f1 = item[0].replace("..",prev_dir)
                if item[1].find("..") != 1:
                    f2 = item[1].replace("..",prev_dir)
                window['-COORDFILE1-'].update(f1)
                window['-COORDFILE2-'].update(f2)
        elif key == 'filetop':
            if isinstance(item, list) and len(item) == 2:
                if item[0].find(".") != 1:
                    f1 = item[0].replace(".",working_dir)
                if item[1].find(".") != 1:
                    f2 = item[1].replace(".",working_dir)
                if item[0].find("..") != 1:
                    f1 = item[0].replace("..",prev_dir)
                if item[1].find("..") != 1:
                    f2 = item[1].replace("..",prev_dir)
                window['-TOPOFILE1-'].update(f1)
                window['-TOPOFILE2-'].update(f2)
        elif key == 'coordination_numbers_Z':
            if item == True:
                window["-ZCheckBox-"].update(True)
                enable_disable_Zoptions(window, dict_Zsamples)
            else:
                window["-ZCheckBox-"].update(False)
                enable_disable_Zoptions(window, dict_Zsamples)
        elif key == 'calculate_volume':
            if item == True:
                window["-VolCheckBox-"].update(True)
            else:
                window["-VolCheckBox-"].update(False)
        elif key == 'interaction_energy':
            if item == True:
                window["-IECheckBox-"].update(True)
                enable_disable_Zoptions(window, dict_QMoptions)
            else:
                window["-IECheckBox-"].update(False)
                enable_disable_Zoptions(window, dict_QMoptions)
        elif key == 'boxl':
            window['-BOXL-'].update(item)
        elif key == 'Z_parameters':
            if isinstance(item, dict):
                for key2, item2 in item.items():
                    if key2 == "Z_samples" and isinstance(item2, int):
                        window['-ZSAMPLES-'].update(item2)
                    elif key2 == "Z_puttrialmonomers" and isinstance(item2, int):
                        window['-ZTRIALS-'].update(item2)
                    elif key2 == 'Z_debug':
                        if item2 == True:
                            window["-ZDEBUG-"].update(True)
                        else:
                            window["-ZDEBUG-"].update(False)
                    elif key2 == "Z_nonbonded" and isinstance(item2, str) and\
                        item2 in window['-ZNONBONDED-'].Values:
                        window['-ZNONBONDED-'].update(item2)
        elif key == 'energy_parameters':
            if isinstance(item, dict):
                for key2, item2 in item.items():
                    if key2 == "qm_engine" and isinstance(item2, str) and \
                            item2 in window['-QMPACK-'].Values:
                        window['-QMPACK-'].update(item2)
                    elif key2 == "qm_path_exe" and isinstance(item2, str):
                        window['-QMPATHEXE-'].update(item2)
                    elif key2 == "qm_scratch_dir" and isinstance(item2, str):
                        window['-QMPATHSCRATCH-'].update(item2)
                    elif key2 == "qm_charge" and isinstance(item2, int):
                        window['-QMCHARGE-'].update(item2)
                    elif key2 == "qm_multiplicity" and isinstance(item2, int):
                        window['-QMMULTIPLICITY-'].update(item2)
                    elif key2 == "qm_basisset" and isinstance(item2, str):
                        window['-QMBASISSET-'].update(item2)
                    elif key2 == "qm_method" and isinstance(item2, str):
                        window['-QMMETHOD-'].update(item2)
                    elif key2 == "qm_task" and isinstance(item2, str):
                        window['-QMTASK-'].update(item2)
                    elif key2 == "qm_otherkeys" and isinstance(item2, str):
                        window['-QMOTHER-'].update(item2)
                    elif key2 == "number_configurations" and isinstance(item2, int):
                        window['-QMCONF-'].update(item2)
                    elif key2 == "qm_memory_mb" and isinstance(item2, int):
                        window['-QMMEM-'].update(item2)
        elif key == 'server':
            if isinstance(item, dict):
                for key2, item2 in item.items():
                    if key2 == "name" and isinstance(item2, str):
                        window['-SERVERNAME-'].update(item2)
                        if item2 == "localhost":
                            window['-RADIO1a-'].update(True)
                            window['-SERVERQUEUE-'].update(disabled=True)
                            window['-SERVERUSER-'].update(disabled=True)
                            window['-SERVERKEYFILE-'].update(disabled=True)
                            window['-BROWSEKEYFILE-'].update(disabled=True)
                            window['-SERVERPARTITION-'].update(disabled=True)
                            window['-SERVERNODELIST-'].update(disabled=True)
                            window['-SERVERPARTMASTER-'].update(disabled=True)
                            window['-SERVERPARTMASTERNODE-'].update(disabled=True)
                    elif key2 == "queue_system" and isinstance(item2, str):
                        window['-SERVERQUEUE-'].update(item2)
                        if item2 == "slurm":
                            window['-RADIO1b-'].update(True)
                            window['-SERVERQUEUE-'].update(disabled=False)
                            window['-SERVERUSER-'].update(disabled=False)
                            window['-SERVERKEYFILE-'].update(disabled=False)
                            window['-BROWSEKEYFILE-'].update(disabled=False)
                            window['-SERVERPARTITION-'].update(disabled=False)
                            window['-SERVERNODELIST-'].update(disabled=False)
                            window['-SERVERPARTMASTER-'].update(disabled=False)
                            window['-SERVERPARTMASTERNODE-'].update(disabled=False)
                    elif key2 == "username" and isinstance(item2, str):
                        window['-SERVERUSER-'].update(item2)
                    elif key2 == "key_file" and isinstance(item2, str):
                        window['-SERVERKEYFILE-'].update(item2)
                    elif key2 == "remote_dir" and isinstance(item2, str):
                        window['-SERVERREMOTEDIR-'].update(item2)
                    elif key2 == "local_dir" and isinstance(item2, str):
                        window['-SERVERLOCALDIR-'].update(item2)
                    elif key2 == "ncpus" and isinstance(item2, int):
                        window['-SERVERCPU-'].update(item2)
                    elif key2 == "partition" and isinstance(item2, str):
                        window['-SERVERPARTITION-'].update(item2)
                    elif key2 == "nodelist" and isinstance(item2, list):
                        l = ""
                        for i in item2:
                            l += i+","
                        l = l[0:-1]
                        window['-SERVERNODELIST-'].update(l)
                    elif key2 == "partitionmaster" and isinstance(item2, str):
                        window['-SERVERPARTMASTER-'].update(item2)
                    elif key2 == "nodelistmaster" and isinstance(item2, str):
                        window['-SERVERPARTMASTERNODE-'].update(item2)
                    elif key2 == "memory" and isinstance(item2, int):
                        window['-SERVERMEMORY-'].update(item2)
        elif key == 'montecarlo_parameters':
            if isinstance(item, dict):
                for key2, item2 in item.items():
                    if key2 == "temperatures":
                        window["-MCINITEMP-"].update(item2[0])
                        window["-MCENDTEMP-"].update(item2[1])
                        window["-MCSTEPTEMP-"].update(item2[2])
                    elif key2 == "mode":
                        window["-MCMODE-"].update(item2)
                    elif key2 == "pairs":
                        window["-MCPAIRS-"].update(item2)
        elif key == 'montecarlo_run':
            if item == True:
                window["-MCCheckBox-"].update(True)
            else:
                window["-MCCheckBox--"].update(False)

    # Delete inputdict
    del(inputdict)

# =============================================================================
def run_chiripa(window):

    script_path = window['-HIDEINPUTSCRIPT-'].get()
    loc_main = window.current_location()
    size = window.size
    loc_center = (loc_main[0]+size[0]/2, loc_main[1]+size[1]/2)
    dir = os.path.dirname(script_path)+"/"
    program = os.path.basename(script_path)

    with open(dir+"tmp.sh", 'w') as f:
        l = "#!/bin/bash\n"
        l += "cd {}\n".format(dir)
        l += "python {}".format(script_path)
        f.writelines(l)

    if os.path.exists(dir+"Z_running_tmp.txt"):
        sg.popup("Z calculation is running", title="Warning", location = loc_center)
        # Check and update database. Do not send anything if exist the database
        proc = subprocess.Popen(['bash', '{}tmp.sh'.format(dir)], close_fds=True)
    elif os.path.exists(dir+"Z_done.txt"):
        sg.popup("Z calculation is already done", title="Warning", location = loc_center)
        # Check and update database. Do not send anything if exist the database
        proc = subprocess.Popen(['bash', '{}tmp.sh'.format(dir)], close_fds=True)
    else:
        # Run Z and QM calculations. The  QM calculations are send only if the database does not exist
        proc = subprocess.Popen(['bash', '{}tmp.sh'.format(dir)], close_fds=True)

# =============================================================================
def run_chiripa_get_qmresults_database(window):

    script_path = window['-HIDEINPUTSCRIPT-'].get()
    dir = os.path.dirname(script_path)+"/"
    program = os.path.basename(script_path)

    fnamedblist = glob.glob(dir + '*.db')
    if len(fnamedblist) != 1:
        return False
    fnamedb = fnamedblist[0]

    # Update the database
    with open(dir+"tmp.sh", 'w') as f:
        l = "#!/bin/bash\n"
        l += "cd {}\n".format(dir)
        l += "python {}".format(script_path)
        f.writelines(l)

    os.system("bash {}tmp.sh".format(dir))

    db  = chi.DBjobs(fnamedb)

    total, inserver, pending, complete, error, running = db.account_qm("qm_jobs")

    Z_val = []; Z_std = []
    if os.path.exists(dir+"Z_results.log"):
        with open(dir+"Z_results.log", "r") as f:
            lines = f.readlines()
            Z_val = ['None',
                     'None',
                     lines[0].split("=")[1].split()[0],
                     lines[1].split("=")[1].split()[0],
                     lines[2].split("=")[1].split()[0],
                     lines[3].split("=")[1].split()[0],]
            Z_std = ['None',
                     'None',
                     lines[0].split("=")[1].split()[2],
                     lines[1].split("=")[1].split()[2],
                     lines[2].split("=")[1].split()[2],
                     lines[3].split("=")[1].split()[2],]

    label1 = "-NAME1-"
    label2 = "-NAME2-"
    name1 = window[label1].get()
    name2 = window[label2].get()
    window["-TABLEQM_0_0-"].update(name1)
    window["-TABLEQM_0_1-"].update('None')
    window["-TABLEQM_0_2-"].update(total[0])
    window["-TABLEQM_0_3-"].update(complete[0])
    window["-TABLEQM_0_4-"].update('None')
    window["-TABLEQM_0_5-"].update('None')
    window["-TABLEQM_1_0-"].update(name2)
    window["-TABLEQM_1_1-"].update('None')
    window["-TABLEQM_1_2-"].update(total[1])
    window["-TABLEQM_1_3-"].update(complete[1])
    window["-TABLEQM_1_4-"].update('None')
    window["-TABLEQM_1_5-"].update('None')

    ind = 2
    for iseg in range(1,3):
        for jseg in range(1,3):
            label1 = "-NAME{0:d}-".format(iseg)
            label2 = "-NAME{0:d}-".format(jseg)
            name1 = window[label1].get()
            name2 = window[label2].get()
            window["-TABLEQM_{}_0-".format(ind)].update(name1)
            window["-TABLEQM_{}_1-".format(ind)].update(name2)
            window["-TABLEQM_{}_2-".format(ind)].update(total[ind])
            window["-TABLEQM_{}_3-".format(ind)].update(complete[ind])
            if len(Z_val) != 0:
                window["-TABLEQM_{}_4-".format(ind)].update(Z_val[ind])
                window["-TABLEQM_{}_5-".format(ind)].update(Z_std[ind])
            ind +=1
    db.close_connection()

# =============================================================================
def run_chiripa_get_mcresults_database(window):

    script_path = window['-HIDEINPUTSCRIPT-'].get()
    dir = os.path.dirname(script_path) + "/"
    program = os.path.basename(script_path)

    foutlist = glob.glob(dir + 'output_results.log')
    if len(foutlist) != 1:
        m = "No MC results to show. Check if 'output_results.log' exists in the working dir.\n"
        window["-MCRESULTS-"].update(m)
        return False
    fout = foutlist[0]

    m = ''
    with open(fout, 'r') as f:
        lines = f.readlines()
        isfound = False
        for iline in lines:
            if iline.find("Monte Carlo results") != -1:
                isfound=True
            if isfound:
                m += iline
                if len(re.sub(r"\s", "", iline)) == 0:
                    break

    window["-MCRESULTS-"].update(m)

    m = ''
    with open(fout, 'r') as f:
        lines = f.readlines()
        isfound = False
        for iline in lines:
            if iline.find("#Flory-Hugging parameters") != -1:
                isfound=True
            if isfound:
                m += iline
                if len(re.sub(r"\s", "", iline)) == 0:
                    break

    window["-CHIRESULTS-"].update(m)

    return True

# =============================================================================
def clean_form(window):

    keys_input_to_clear = ['-NAME1-', '-NAME2-', '-COORDFILE1-', '-COORDFILE2-', '-TOPOFILE1-', '-TOPOFILE2-',
                           '-ZSAMPLES-', '-ZTRIALS-','-QMPATHEXE-', '-QMPATHSCRATCH-', '-QMMEM-', '-QMCONF-',
                           '-QMCHARGE-', '-QMMULTIPLICITY-', '-QMMETHOD-', '-QMBASISSET-', '-QMOTHER-',
                           '-SERVERNAME-', '-SERVERQUEUE-', '-SERVERUSER-', '-SERVERKEYFILE-', '-SERVERREMOTEDIR-',
                           '-SERVERLOCALDIR-', '-SERVERCPU-', '-SERVERPARTITION-', '-SERVERNODELIST-',
                           '-SERVERPARTMASTER-', '-SERVERPARTMASTERNODE-', '-SERVERMEMORY-']

    keys_checkbox_to_clear = {"-TOPOCHECKBOX1-": True, "-TOPOCHECKBOX2-": True,
                              "-ZCheckBox-": False, "-VolCheckBox-": False, "-IECheckBox-": False,
                              '-ZDEBUG-': False}

    keys_combo_to_clear = {'-ZNONBONDED-': 'truhlar', '-QMPACK-':'gaussian','-QMTASK-':'energy'}

    for ikey in keys_input_to_clear:
        window[ikey].update("")
        for irow in range(6):
            for icol in range(6):
                window["-TABLEQM_{}_{}-".format(irow, icol)].update("")

    for ikey, item in keys_checkbox_to_clear.items():
        window[ikey].update(item)
        enable_disable_topology(window, dict_topology)

    for ikey, item in keys_combo_to_clear.items():
        window[ikey].update(item)

    window['-RADIO1a-'].update(True)
    window['-SERVERQUEUE-'].update("", disabled=True)
    window['-SERVERUSER-'].update("", disabled=True)
    window['-SERVERKEYFILE-'].update("", disabled=True)
    window['-BROWSEKEYFILE-'].update(disabled=True)
    window['-SERVERPARTITION-'].update("", disabled=True)
    window['-SERVERNODELIST-'].update("", disabled=True)
    window['-SERVERPARTMASTER-'].update("", disabled=True)
    window['-SERVERPARTMASTERNODE-'].update("", disabled=True)

# =============================================================================
def waiting_for_events(window, event, values):
    """
    All events in the GUI
    """

    filename_import = None

    # ================== BUTTONVIEWMOL ==================
    if event == '-BUTTONVIEWMOL-':
        try:
            moleculepath = [values['-COORDFILE1-'], values['-COORDFILE2-']]
            molecule_web(window, moleculepath)
        except:
            print (moleculepath)
            s = window.size
            l = window.current_location()
            sg.popup("ERROR!!!!", "Coordinates file of Segment 1 or 2 either does not exist or the format is incorrect.",
                     location=(l[0]+s[0]/2.0, l[1]+s[1]/2.0))

    if event == '-HIDEBUTTON1-' or event == "Tab:23":
        enable_disable_buttonscript(window, dict_button_script, dict_Zsamples, dict_Serveroptions)

    # ================== TOPOLOGY OBJECTS ==================
    if event == "-TOPOCHECKBOX1-" or event == "-TOPOCHECKBOX2-":
        enable_disable_topology(window, dict_topology)

    # ================== COORDINATE FILE OBJECTS ==================
    for key, items in dict_coord.items():
        try:
            newvalue = values[key]
            if newvalue != '' and event == key :
                name = newvalue.split("/")[-1].split(".")[0]
                window[items[1]].update(name)
                if window[items[2]].get():
                    window[items[3]].update(newvalue)
        except:
            pass

    # ================== Z Options ==================
    if event == '-ZCheckBox-':
        enable_disable_Zoptions(window, dict_Zsamples)

    # ================== TOPOLOGY FILE OBJECTS ==================
    #####enable_disable_buttonscript(window, dict_button_script, dict_Zsamples, dict_Serveroptions)
    if event == '-BUTTONCREATESCRIPT-':
        loc = window.CurrentLocation()
        window.hide()
        l = create_script_chiripa(window)
        winaux = sg.Window("Python script", [[sg.Multiline(l,size=(150,40))]], location=loc )
        while True:
            event1, values1 = winaux.read()
            if event1 == sg.WIN_CLOSED:
                window.un_hide()
                break

    # ============= SERVER LAYOUT EVENTS =============
    if event == '-RADIO1a-':
        window['-SERVERNAME-'].update(value="localhost")
        window['-SERVERQUEUE-'].update("", disabled=True)
        window['-SERVERUSER-'].update("", disabled=True)
        window['-SERVERKEYFILE-'].update("", disabled=True)
        window['-BROWSEKEYFILE-'].update(disabled=True)
        window['-SERVERPARTITION-'].update("",disabled=True)
        window['-SERVERNODELIST-'].update("", disabled=True)
        window['-SERVERPARTMASTER-'].update("", disabled=True)
        window['-SERVERPARTMASTERNODE-'].update("", disabled=True)
        window['-SERVERMEMORY-'].update("", disabled=False)

    if event == '-RADIO1b-':
        window['-SERVERNAME-'].update(value="trueno.csic.es")
        window['-SERVERQUEUE-'].update("slurm", disabled=False)
        window['-SERVERUSER-'].update("jramos", disabled=False)
        window['-SERVERKEYFILE-'].update("/home/jramos/.ssh/id_rsa_chiripa", disabled=False)
        window['-BROWSEKEYFILE-'].update(disabled=False)
        window['-SERVERPARTITION-'].update("simacro", disabled=False)
        window['-SERVERNODELIST-'].update("trueno90, trueno91, trueno92", disabled=False)
        window['-SERVERPARTMASTER-'].update("simacro", disabled=False)
        window['-SERVERPARTMASTERNODE-'].update("trueno36", disabled=False)
        window['-SERVERMEMORY-'].update("", disabled=False)

    # ============= Menu Events =============
    if event == "Import...":
        clean_form(window)
        loc = window.current_location()
        filename = sg.popup_get_file('Import file to GUI', no_window=False,
                                     location=loc, initial_folder="../tutorials/")

        import_file_to_gui(window, filename)
        enable_disable_buttonscript(window, dict_button_script, dict_Zsamples, dict_Serveroptions)
        window['-HIDEINPUTSCRIPT-'].update(filename)

    if event == "Export...":
        loc = window.current_location()
        filename = sg.popup_get_file('Export file from GUI', no_window=False,
                                     location=loc, initial_folder="./", save_as=True)
        create_script_chiripa(window, save=True, filename=filename)
        window['-HIDEINPUTSCRIPT-'].update(filename)

    if event == "Clean Form":
        clean_form(window)
        window['-HIDEINPUTSCRIPT-'].update('')

    if  event == "About...":
        loc = window.current_location()
        window.disappear()
        sg.popup('About this program', 'Version 1.0',
        'PySimpleGUI Version', sg.version, grab_anywhere=True, location=loc)
        window.reappear()

    # ============= BUTTONRUN: Run Chiripa =============
    if event == '-BUTTONRUN-':
        script_path = window['-HIDEINPUTSCRIPT-'].get()
        dir = os.path.dirname(script_path) + "/"
        fnamedblist = glob.glob(dir + '*.db')
        if len(fnamedblist) == 0:
            run_chiripa(window)
        else:
            loc = window.current_location()
            sg.popup("QM calculations seems to be sent. Try GetResults Button",  grab_anywhere=True, location=loc)

    # ============= BUTTONGETRESULTS: Get Chiripa Results =============
    if event == '-BUTTONGETRESULTS-':
        run_chiripa_get_qmresults_database(window)
        run_chiripa_get_mcresults_database(window)

    # ============= Summary of MC calculations: Save data =============
    if event == "-MCSAVE2-":
        fpathname = values["-MCSAVE2-"]
        with open(fpathname, 'w') as f:
            m = values['-MCRESULTS-']
            f.write(m)
            mm = "File {} has been saved!!!".format(fpathname)
            window['-MCTEXT-'].update(value=mm, visible=True)

    # ============= Summary of CHI calculations: Save data =============
    if event == "-CHISAVE2-":
        fpathname = values["-CHISAVE2-"]
        with open(fpathname, 'w') as f:
            m = values['-CHIRESULTS-']
            f.write(m)
            mm = "File {} has been saved!!!".format(fpathname)
            window['-CHITEXT-'].update(value=mm, visible=True)
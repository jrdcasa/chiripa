import PySimpleGUI as sg
from common_elements import output_box ,set_options, theme

msg1 = "Path to the QM executable in the server."
msg2 = "Path to the QM scratcch dir in the server."
msg3 = "Memory in Mbytes to run the QM calculation."
msg4 = "Number of different configrurations for segment pairs"

sg.ChangeLookAndFeel(theme)

col1 = sg.Column([
        [sg.Frame('QM package details',[
            [sg.Text("qm_engine : ", size=(14, 1), pad=((0,0),(20,0))),
             sg.Combo(['gaussian', 'gamess', 'nwchem'], disabled=False, key='-QMPACK-',
                      size=(40, 1), pad=((0,0),(20,0)))],
            [sg.Text("qm_path_exe : ", size=(14, 1), pad=((0,0),(20,0))),
             sg.Input(key="-QMPATHEXE-", tooltip=msg1, disabled=False, size=(50, 1), pad=((0,0),(20,0)))],
            [sg.Text("qm_scratch_dir : ", size=(14, 1), pad=((0,0),(20,30))),
             sg.Input(key="-QMPATHSCRATCH-", tooltip=msg2, disabled=False, size=(50, 1), pad=((0,0),(20,30)))]
        ], pad=((30,10),(10,0)), title_color='blue'),
        ]])

col2 = sg.Column([
        [sg.Frame('QM options',[
            [sg.Text("qm_memory_mb : ", size=(25, 1), pad=((0,0),(20,0))),
             sg.Input(key="-QMMEM-", disabled=True, tooltip=msg3, size=(10, 1), pad=((0,0),(20,0)))],
            [sg.Text("qm_number_configurations : ", size=(25, 1), pad=((0,0),(20,20))),
             sg.Input(key="-QMCONF-", disabled=False, tooltip=msg4, size=(10, 1), pad=((0,0),(20,20)))]
        ], pad=((70,10),(20,0)),title_color='blue'),
        ]])

frame1 = sg.Frame('Definition of the QM calculation', [
            [sg.Text('qm_charge : ', size=(12, 1), pad=((10,10),(20,0))),
             sg.Input(key='-QMCHARGE-', size=(5,1), disabled=False, pad=((2,4),(20,0))),
             sg.Text('qm_multiplicity : ', size=(14, 1), pad=((20,0),(20,0))),
             sg.Input(key='-QMMULTIPLICITY-', size=(5,1), disabled=False, pad=((2,4),(20,0))),
             sg.Text('qm_Task : ', size=(12, 1), pad=((20,0),(20,0))),
             sg.Combo(['energy', 'optimization'], disabled=False, key='-QMTASK-', pad=((2, 4),(20,0))) ],
            [sg.Text('qm_method : ', size=(12, 1), pad=((10,0),(20,0))),
             sg.Input(key='-QMMETHOD-', size=(12,1), disabled=False, pad=((2,4),(20,0))),
             sg.Text('qm_basisset : ', size=(12, 1), pad=((120,0),(20,0))),
             sg.Input(key='-QMBASISSET-', size=(12,1), disabled=False, pad=((2, 4),(20,0))) ],
            [sg.Text('qm_otherkeys : ', size=(12, 1), pad=((10,10),(20,20))),
             sg.Input(key='-QMOTHER-', size=(50, 1), disabled=False, pad=((2,4),(20,20))),],
         ], pad=((30,10), (10,0)),title_color='blue')

qm_layout       = [[col1, col2], [frame1]]




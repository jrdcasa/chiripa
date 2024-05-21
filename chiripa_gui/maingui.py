import PySimpleGUI as sg
import os
from main_tabgroup import tabgroup_layout
from menu_layout import menu_layout
from events import waiting_for_events
from common_elements import output_box, row_buttons, theme

# MAIN WINDOW LOOP
sg.ChangeLookAndFeel(theme)
layout =[ [menu_layout], [ tabgroup_layout ], [row_buttons]]
window = sg.Window('chiripa GUI', layout, size = (1300,750), location=(100,100),
                   finalize=True, return_keyboard_events=True, force_toplevel=False)
window.Element("-TABGROUP-").set_size(size=(1300, 600))
window.Element("-TABGROUP_RES-").set_size(size=(1300, 600))

while True:

    event, values = window.read()

    if event == sg.WIN_CLOSED or event == "Exit":
        script_path = window['-HIDEINPUTSCRIPT-'].get()
        dir = os.path.dirname(script_path)+"/"
        if os.path.exists(dir+"tmp.sh"):
            os.remove(dir+"tmp.sh")
        break

    waiting_for_events(window, event, values)
    #print(event, values)








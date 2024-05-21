import PySimpleGUI as sg
from common_elements import set_options, theme

sg.ChangeLookAndFeel(theme)

# ------ Menu Definition ------ #
menu_def = [['File', ['Import...', '!Export...', 'Exit']],
            ['Edit', ['Clean Form', ], ],
            ['Help', 'About...'], ]

menu_layout = sg.Menu(menu_def, tearoff=False, pad=(200, 1),
                      font=("Helvetica", 14), key="-MENU-")
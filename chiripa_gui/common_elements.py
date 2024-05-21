import PySimpleGUI as sg

theme="GreenMono"
set_options = sg.set_options(text_justification='left', font=("Helvetica", 13))
sg.ChangeLookAndFeel(theme)

# OUTPUT BOX is common to all folders in the GUI
output_box = sg.Frame(title = "Output", layout=[
        [sg.Multiline(size=(140,10), key='-OUTPUTBOX-', autoscroll=True)]
       ], title_color='blue', pad=((0,10), (30,0)))

row_buttons = sg.Column ([
        [sg.Button('View Segments', key='-BUTTONVIEWMOL-'),
         sg.Button('Run Chiripa', key='-BUTTONRUN-', disabled=True),
         sg.Button('Create Python Script', key='-BUTTONCREATESCRIPT-', disabled=True),
         sg.Button('Get Results', key='-BUTTONGETRESULTS-', disabled=False)]
   ], pad=((650,10), (30,00)))



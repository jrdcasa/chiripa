import PySimpleGUI as sg
from common_elements import set_options


col1 = sg.Column([
        [sg.Image("/home/jramos/Imagenes/ZCoordination_small.png", pad=((180, 10),(20, 0)))]
     ])

col2 = sg.Column([
        [sg.Text('Z samples      :', size=(14,1), pad=((0,0),(20,0)), font=("Helvetica", 14)),
         sg.Input(key='-ZSAMPLES-', size=(5,1), disabled=True, pad=((20,0),(20,0)), font=("Helvetica", 14),
                  enable_events=True),
         sg.Button('OK', key="-HIDEBUTTON1-", bind_return_key=True, visible=False, pad=((20,0),(20,0)))],
        [sg.Text('Z trials       :', size=(14,1), pad=((0,0),(20, 0)), font=("Helvetica", 14)),
         sg.Input(key='-ZTRIALS-', size=(5,1), disabled=True, pad=((20,0),(20,0)), font=("Helvetica", 14),
                  enable_events=True )],
        [sg.Text('Z debug        :', size=(14,1), pad=((0,0),(10,0)), font=("Helvetica", 14)),
         sg.Checkbox('', key='-ZDEBUG-', size=(15,1), disabled=True, default=False, pad=((20,0),(20,0)),
                     font=("Helvetica", 14), enable_events=True)],
        [sg.Text('Z non-bonded   :', size=(14, 1), pad=((0,0),(20,0)), font=("Helvetica", 14)),
         sg.Combo(['truhlar', 'okuwaki'], readonly=True, disabled=True, default_value='truhlar',
                  key='-ZNONBONDED-',size=(15, 1), pad=((20,0),(20,0)), font=("Helvetica", 14), enable_events=True)],
    ], pad=((200,0),(10,0)))

zoptions_layout       = [[col1, col2]]

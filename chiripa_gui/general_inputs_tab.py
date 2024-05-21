import PySimpleGUI as sg
from common_elements import output_box, set_options, theme

sg.ChangeLookAndFeel(theme)

col1 = sg.Column([
    # Information frame
    [sg.Frame('Segment 1:', [[sg.Text(),
                              sg.Column([[sg.Text('Name', size=(14, 1)),
                                          sg.Input(key='-NAME1-', size=(26, 1),
                                                   tooltip='Name of Segment 1')],
                                         [sg.Text('Coordinate File:', size=(14, 1)),
                                          sg.Input(key='-COORDFILE1-', size=(26, 1),
                                                   tooltip='Coordinate file of Segment 1', enable_events=True),
                                          sg.FileBrowse(button_text="Browse", key="-BROWSECOORD1-")],
                                         [sg.Text('Topology File:', size=(14, 1)),
                                          sg.Input(key='-TOPOFILE1-', size=(26, 1),
                                                   tooltip='Topology file of Segment 1', disabled=True,
                                                   enable_events=True),
                                          sg.FileBrowse(button_text="Browse", key="-BROWSETOPO1-", disabled=True), ],
                                         [sg.Checkbox("Topology?", key="-TOPOCHECKBOX1-", default=True,
                                                      enable_events=True)],
                                         ]
                                        )]], title_color='blue', pad=((30, 10), (30, 0))
              )]
])

col2 = sg.Column([
    # Information frame
    [sg.Frame('Segment 2:', [[sg.Text(),
                              sg.Column([[sg.Text('Name:', size=(14, 1)),
                                          sg.Input(key='-NAME2-', size=(26, 1), tooltip='Name of Segment 2')],
                                         [sg.Text('Coordinate File:', size=(14, 1)),
                                          sg.Input(key='-COORDFILE2-', size=(26, 1),
                                                   tooltip='Coordinate file of Segment 2', enable_events=True),
                                          sg.FileBrowse(button_text="Browse", key="-BROWSECOORD2-")],
                                         [sg.Text('Topology File:', size=(14, 1)),
                                          sg.Input(key='-TOPOFILE2-', size=(26, 1),
                                                   tooltip='Topology file of Segment 2',
                                                   disabled=True, enable_events=True),
                                          sg.FileBrowse(button_text="Browse", key="-BROWSETOPO2-", disabled=True)],
                                         [sg.Checkbox("Topology?", key="-TOPOCHECKBOX2-",
                                                      default=True, enable_events=True)]
                                         ],
                                        )]], title_color='blue', pad=((10, 10), (30, 0))
              )]
])

checkboxes_inputs = sg.Frame(title="Calculations", layout=[
    [sg.Checkbox('Coordination Number (Z-parameter)', key="-ZCheckBox-", size=(45, 2),
                 pad=((0, 0), (10, 10)), enable_events=True, font=("Helvetica", 16)),
     sg.Checkbox('Calculate Volume', key="-VolCheckBox-", size=(37, 2),
                 pad=((0, 0), (10, 10)), enable_events=True, font=("Helvetica", 16))],
    [sg.Checkbox('Calculate Interaction Energies', key="-IECheckBox-", size=(45, 2),
                 pad=((0, 0), (10, 10)), enable_events=True, font=("Helvetica", 16)),
     sg.Checkbox('Monte Carlo Run', key="-MCCheckBox-", size=(37, 2),
                 pad=((0, 0), (10, 10)), enable_events=True, font=("Helvetica", 16))
     ],
], title_color='blue', pad=((30, 10), (10, 0)))

box_inputs = sg.Frame(title="Box", layout=[
    [sg.Text('Box length (A) :', size=(14, 1), pad=((10, 10), (20, 20)), font=("Helvetica", 14)),
     sg.Input(key='-BOXL-', size=(5, 1), disabled=False, pad=((10, 900), (20, 20)), font=("Helvetica", 14),
              enable_events=True, default_text=30)],
], title_color='blue', pad=((30, 10), (10, 0)))

hide_input = sg.Input(key='-HIDEINPUTSCRIPT-', visible=False)

# Change font and text of all TEXT ELEMENTS. It seems that the FONT defined in set_options does not take effect
general_inputs_layout = [[col1, col2], [checkboxes_inputs], [box_inputs], [hide_input]]

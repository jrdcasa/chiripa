import PySimpleGUI as sg

"""
    Monte Carlo Tab
"""

col = sg.Column([
    [sg.Text('Initial temperature (K) :', size=(24, 1), pad=((0, 0), (20, 0)), font=("Helvetica", 14)),
     sg.Input(key='-MCINITEMP-', size=(5, 1), pad=((20, 0), (20, 0)), font=("Helvetica", 14),
              enable_events=True, default_text=300)],
    [sg.Text('End Temperature (K) :', size=(24, 1), pad=((0, 0), (20, 0)), font=("Helvetica", 14)),
     sg.Input(key='-MCENDTEMP-', size=(5, 1), pad=((20, 0), (20, 0)), font=("Helvetica", 14),
              enable_events=True, default_text=300)],
    [sg.Text('Step Temperature (K) :', size=(24, 1), pad=((0, 0), (20, 0)), font=("Helvetica", 14)),
     sg.Input(key='-MCSTEPTEMP-', size=(5, 1), pad=((20, 0), (20, 0)), font=("Helvetica", 14),
              enable_events=True, default_text=1)],
    [sg.Text('Monte Carlo Mode :', size=(24, 1), pad=((0, 0), (20, 0)), font=("Helvetica", 14)),
     sg.Combo(['sequential'], readonly=True, default_value='sequential',
              key='-MCMODE-', size=(15, 1), pad=((20, 0), (20, 0)), font=("Helvetica", 14),
              enable_events=True)],
    [sg.Text('Pairs MC :', size=(24, 1), pad=((0, 0), (20, 0)), font=("Helvetica", 14)),
     sg.Combo(['all'], readonly=True, default_value='all',
              key='-MCPAIRS-', size=(15, 1), pad=((20, 0), (20, 0)), font=("Helvetica", 14),
              enable_events=True)],
    ], pad=((200, 0), (10, 0)))

mc_layout = [[col]]

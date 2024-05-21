import PySimpleGUI as sg


headings = ['Segment 1', 'Segment 2', "Total QM send", "Completed QM jobs", "Z-value", "Std Z-Value"]

frame1 = sg.Column([
    [sg.Frame('Summary of QM calculations', [
        [sg.Text(' ')] + [sg.Text(h, size=(16,1),  pad=(3,0), justification="center") for h in headings],
        [sg.Input(size=(16,1), disabled=True, justification="center", key="-TABLEQM_{}_{}-".format(0,col)) for col in range(6)],
        [sg.Input(size=(16,1), disabled=True, justification="center", key="-TABLEQM_{}_{}-".format(1,col)) for col in range(6)],
        [sg.Input(size=(16,1), disabled=True, justification="center", key="-TABLEQM_{}_{}-".format(2,col)) for col in range(6)],
        [sg.Input(size=(16,1), disabled=True, justification="center", key="-TABLEQM_{}_{}-".format(3,col)) for col in range(6)],
        [sg.Input(size=(16,1), disabled=True, justification="center", key="-TABLEQM_{}_{}-".format(4,col)) for col in range(6)],
        [sg.Input(size=(16,1), disabled=True, justification="center", key="-TABLEQM_{}_{}-".format(5,col)) for col in range(6)],
        [sg.Text(' ')],]
    , pad=((20,0),(20,20)), title_color='blue')]
    ])

qm_results_layout = [[frame1]]
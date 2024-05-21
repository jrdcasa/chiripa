import PySimpleGUI as sg


frame1 = sg.Column([
    [sg.Frame('Summary of MC calculations', [
         [sg.Multiline(size=(88, 20), disabled=True, key="-MCRESULTS-")],
    ], pad=((20, 0), (20, 20)), title_color='blue')]
     ])

frame2_col1 = sg.Column([
    [sg.Frame('', [
        [sg.InputText(visible=False, enable_events=True, key='-MCSAVE1-')],
        [sg.FileSaveAs(key='-MCSAVE2-',
            enable_events=True,
            file_types=(("Text Files", "*.txt"),),
        )],
    ], pad=((20, 0), (20, 20)), title_color='blue')]
])

frame2_col2 = sg.Column([

         [sg.Text(text = "", visible=False, key='-MCTEXT-', size=(90, 10), pad=((60, 0), (40, 20)) )],


 ])

mc_results_layout = [[frame1], [frame2_col1, frame2_col2]]
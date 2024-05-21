import PySimpleGUI as sg

frame1 = sg.Column([
    [sg.Frame('Summary of CHI parameters', [
         [sg.Multiline(size=(88, 20), disabled=True, key="-CHIRESULTS-")],
    ], pad=((20, 0), (20, 20)), title_color='blue')]
     ])

frame2_col1 = sg.Column([
    [sg.Frame('', [
        [sg.InputText(visible=False, enable_events=True, key='-CHISAVE1-')],
        [sg.FileSaveAs(key='-CHISAVE2-',
            enable_events=True,
            file_types=(("Text Files", "*.txt"),),
        )],
    ], pad=((20, 0), (20, 20)), title_color='blue')]
])

frame2_col2 = sg.Column([

         [sg.Text(text = "", visible=False, key='-CHITEXT-', size=(90, 10), pad=((60, 0), (40, 20)) )],


 ])


chi_results_layout = [[frame1], [frame2_col1, frame2_col2]]
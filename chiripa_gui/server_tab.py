import PySimpleGUI as sg
from common_elements import theme


sg.ChangeLookAndFeel(theme)

choices = ["localhost", "slurm"]

col1 = sg.Column([
            [sg.Frame('Choose server type', [
                [sg.Radio(choices[0], "-RADIO1-", key="-RADIO1a-", default=True, enable_events=True, pad=((0,10)))],
                [sg.Radio(choices[1], "-RADIO1-", key="-RADIO1b-", enable_events=True, pad=((0,0),(0,10)))],
        ], title_color='blue')]], pad=((30,0), (125,0)), size=(175, 125))

col2 = sg.Column([
            [sg.Frame('', [
                [sg.Text("Server Name :", size=(12,1), pad=(10,10)),
                 sg.Input(key="-SERVERNAME-", size=(20,1), default_text="localhost", pad=(10,10))],
                [sg.Text("Queue system :", size=(12,1), pad=(10,10)),
                 sg.Input(key="-SERVERQUEUE-", size=(20,1), default_text="", disabled=True, pad=(10,10))],
                [sg.Text("Username :", size=(12, 1), pad=(10,10)),
                 sg.Input(key="-SERVERUSER-", size=(20, 1), default_text="", disabled=True, pad=(10,10))],
                [sg.Text("Key file :", size=(12, 1), pad=(10,10)),
                 sg.Input(key="-SERVERKEYFILE-", size=(20, 1), default_text="", disabled=True, pad=(10,10)),
                 sg.FileBrowse(button_text = "Browse", key="-BROWSEKEYFILE-", disabled=True, pad=(10,10))],
                [sg.Text("Remote dir :", size=(12, 1), pad=(10,10)),
                 sg.Input(key="-SERVERREMOTEDIR-", size=(20, 1), default_text="", pad=(10,10)),
                 sg.FileBrowse(button_text="Browse", key="-BROWSEREMOTEDIR-", pad=(10,10))],
                [sg.Text("Local dir :", size=(12, 1), pad=(10,10)),
                 sg.Input(key="-SERVERLOCALDIR-", size=(20, 1), default_text="", pad=(10,10)),
                 sg.FileBrowse(button_text="Browse", key="-BROWSELOCALDIR-", pad=(10,10))],
            ])]], pad=((5,0), (20,00)), size=(500, 300))

col3 = sg.Column([
            [sg.Frame('', [
                [sg.Text("Number of cpus :", size=(18,1), pad=(10,10)),
                 sg.Input(key="-SERVERCPU-", size=(20,1), default_text="1", pad=(10,10))],
                [sg.Text("Partition :", size=(18,1), pad=(10,10)),
                 sg.Input(key="-SERVERPARTITION-", size=(20,1), default_text="", disabled=True, pad=(10,10))],
                [sg.Text("Nodelist :", size=(18, 1), pad=(10,10)),
                 sg.Input(key="-SERVERNODELIST-", size=(20, 1), default_text="", disabled=True, pad=(10,10))],
                [sg.Text("Partition master :", size=(18, 1), pad=(10,10)),
                 sg.Input(key="-SERVERPARTMASTER-", size=(20, 1), default_text="", disabled=True, pad=(10,10)),],
                [sg.Text("Nodelist master :", size=(18, 1), pad=(10,10)),
                 sg.Input(key="-SERVERPARTMASTERNODE-", size=(20, 1), default_text="", disabled=True, pad=(10,10)),],
                [sg.Text("Memory per task (Mb):", size=(18, 1), pad=(10, 10)),
                 sg.Input(key="-SERVERMEMORY-", size=(20, 1), default_text="", disabled=False, pad=(10, 10)), ],
                [sg.Text("Maximum Queue Jobs:", size=(18, 1), pad=(10, 10)),
                 sg.Input(key="-SERVERMAXQUEUEJOBS-", size=(20, 1), default_text="40", disabled=False, pad=(10, 10)), ],
            ])]], pad=((5,0), (70,00)), size=(500, 350))

server_layout   = [[col1, col2 , col3]]

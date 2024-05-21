import PySimpleGUI as sg
from general_inputs_tab import general_inputs_layout
from zoptions_tab import zoptions_layout
from qmlayout_tab import qm_layout
from common_elements import theme
from server_tab import server_layout
from result_tab import result_layout
from mc_tab import mc_layout

sg.ChangeLookAndFeel(theme)

tabgroup_layout = sg.TabGroup([[
                  sg.Tab('General Inputs', general_inputs_layout),
                  sg.Tab('Z-parameter Options', zoptions_layout),
                  sg.Tab('Energy Calc. Options', qm_layout),
                  sg.Tab('Server Options', server_layout),
                  sg.Tab('Monte Carlo parameters', mc_layout),
                  sg.Tab('Results', result_layout),
                ]], key="-TABGROUP-")


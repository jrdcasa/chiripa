import PySimpleGUI as sg
from common_elements import theme
from qm_results_tab import qm_results_layout
from mc_results_tab import mc_results_layout
from chi_results_tab import chi_results_layout

sg.ChangeLookAndFeel(theme)

result_tab_layout = sg.TabGroup([[
                     sg.Tab('QM calculations Summary',qm_results_layout),
                     sg.Tab('Monte Carlo Summary', mc_results_layout),
                     sg.Tab('Chi Parameters', chi_results_layout),
]], key='-TABGROUP_RES-')



result_layout = [[result_tab_layout]]

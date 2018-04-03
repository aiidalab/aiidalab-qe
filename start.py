import ipywidgets as ipw

def get_start_widget(appbase, jupbase, notebase):
    #http://fontawesome.io/icons/
    template = """
    <table>
    <tr>
        <th style="text-align:center">Quantum Espresso</th>
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Cp2k</th>
    <tr>
    <td valign="top"><ul>
    <li><a href="{appbase}/qe/cell_opt.ipynb" target="_blank">Optimize unitcell</a>
    <li><a href="{appbase}/qe/equation_of_states.ipynb" target="_blank">Equation of states</a>
    </ul></td>
    
    <td valign="top"><ul>
    <li><a href="{appbase}/cp2k/cell_opt.ipynb" target="_blank">Optimize unitcell</a>
    <li><a href="{appbase}/cp2k/equation_of_states.ipynb" target="_blank">Equation of states</a>
    </ul></td>

    </tr></table>
"""
    
    html = template.format(appbase=appbase, jupbase=jupbase, notebase=notebase)
    return ipw.HTML(html)
    
#EOF

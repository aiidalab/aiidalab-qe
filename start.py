import ipywidgets as ipw

def get_start_widget(appbase, jupbase):
    #http://fontawesome.io/icons/
    template = """
    <table>
    <tr>
        <th style="text-align:center">Quantum Espresso</th>
        <th style="width:70px" rowspan=2></th>
        <th style="text-align:center">Cp2k</th>        
    <tr>
    <td valign="top"><ul>
    <li><a href="{appbase}/upload_structure.ipynb" target="_blank">Optimize unitcell</a>
    </ul></td>
    
    <td valign="top"><ul>
    <li><a href="{appbase}/nanoribbon/submit.ipynb" target="_blank">Optimize unitcell</a>
    </ul></td>
    </tr></table>
"""
    
    html = template.format(appbase=appbase, jupbase=jupbase)
    return ipw.HTML(html)
    
#EOF

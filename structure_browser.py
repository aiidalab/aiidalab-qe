from aiida.orm.querybuilder import QueryBuilder
from aiida.orm.data.structure import StructureData
from aiida.orm.data.cif import CifData
from aiida.orm.calculation.job import JobCalculation
from aiida.orm.calculation.work import WorkCalculation
from aiida.orm import Node
    
from collections import OrderedDict
import ipywidgets as ipw
import datetime

class StructureBrowser(ipw.VBox):
    
    def __init__(self):
        # Find all process labels
        qb = QueryBuilder()
        qb.append(WorkCalculation,
                  project="attributes._process_label",
                  filters={
                      'attributes': {'!has_key': 'source_code'}
                  }
        )
        qb.order_by({WorkCalculation:{'ctime':'desc'}})
        process_labels = []
        for i in qb.iterall():
            if i[0] not in process_labels:
                process_labels.append(i[0])

        layout = ipw.Layout(width="900px")

        self.mode = ipw.RadioButtons(options=['all', 'uploaded', 'edited', 'calculated'],
                                     layout=ipw.Layout(width="25%"))
        
        
        # Date range
        self.dt_now = datetime.datetime.now()
        self.dt_end = self.dt_now - datetime.timedelta(days=7)
        self.date_start = ipw.Text(value='',
                                   description='From: ',
                                   style={'description_width': '120px'})

        self.date_end = ipw.Text(value='',
                                 description='To: ')

        self.date_text = ipw.HTML(value='<p>Select the date range:</p>')
        
        self.btn_date = ipw.Button(description='Search',
                                   layout={'margin': '1em 0 0 0'})

        self.age_selection = ipw.VBox([self.date_text, ipw.HBox([self.date_start, self.date_end]), self.btn_date],
                                      layout={'border': '1px solid #fafafa', 'padding': '1em'})



        self.btn_date.on_click(self.search)
        self.mode.observe(self.search, names='value')
        
        hr = ipw.HTML('<hr>')
        box = ipw.VBox([self.age_selection,
                        hr,
                        ipw.HBox([self.mode])])
        
        self.results = ipw.Dropdown(layout=layout)
        self.search()
        super(StructureBrowser, self).__init__([box, hr, self.results])
    
    
    def preprocess(self):
        qb = QueryBuilder()
        filters = {'extras': {'!has_key': 'formula'}}
        filters['or'] = [{'type':CifData._plugin_type_string},{'type':StructureData._plugin_type_string}]
        qb.append(Node, filters=filters)
        for n in qb.all(): # iterall() would interfere with set_extra()
            try:
                formula = n[0].get_formula()
            except:
                formula = n[0].get_formulae()[0]
            n[0].set_extra("formula", formula)

    
    def search(self, c=None):
        self.preprocess()
        
        qb = QueryBuilder()
        try: # If the date range is valid, use it for the search
            self.start_date = datetime.datetime.strptime(self.date_start.value, '%Y-%m-%d')
            self.end_date = datetime.datetime.strptime(self.date_end.value, '%Y-%m-%d') + datetime.timedelta(hours=24)
        except ValueError: # Otherwise revert to the standard (i.e. last 7 days)
            self.start_date = self.dt_end
            self.end_date = self.dt_now + datetime.timedelta(hours=24)

            self.date_start.value = self.start_date.strftime('%Y-%m-%d')
            self.date_end.value = self.end_date.strftime('%Y-%m-%d')

        filters = {}
        filters['ctime'] = {'and':[{'<=': self.end_date},{'>': self.start_date}]}
        filters['or'] = [{'type':CifData._plugin_type_string},{'type':StructureData._plugin_type_string}]
        if self.mode.value == "uploaded":
            qb2 = QueryBuilder()
            qb2.append(StructureData, project=["id"])
            qb2.append(Node, input_of=StructureData)
            processed_nodes = [n[0] for n in qb2.all()]
            if processed_nodes:
                filters['id'] = {"!in":processed_nodes}
            qb.append(StructureData, filters=filters)

        elif self.mode.value == "calculated":
            qb.append(JobCalculation)
            qb.append(Node, output_of=JobCalculation, filters=filters)

        elif self.mode.value == "edited":
            qb.append(WorkCalculation)
            qb.append(StructureData, output_of=WorkCalculation, filters=filters)

        else:
            self.mode.value == "all"
            qb.append(Node, filters=filters)

#        qb.order_by({StructureData:{'ctime':'desc'}})
        matches = set([n[0] for n in qb.iterall()])
        matches = sorted(matches, reverse=True, key=lambda n: n.ctime)
        
        c = len(matches)
        options = OrderedDict()
        options["Select a Structure (%d found)"%c] = False

        for n in matches:
            label  = "PK: %d" % n.pk
            label += " | " + n.ctime.strftime("%Y-%m-%d %H:%M")
            label += " | " + n.get_extra("formula")
            label += " | " + n.description
            options[label] = n

        self.results.options = options

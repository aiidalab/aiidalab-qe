import json
from urllib.parse import urlsplit, parse_qs

import ase
import requests
import traitlets
import ipywidgets as ipw

from osp.core.utils import deserialize


class CUDSImporter(ipw.HBox):
    
    title = 'CUDS'
    structure = traitlets.Instance(ase.Atoms, allow_none=True)
    
    def __init__(self, **kwargs):
        self.file_upload = ipw.FileUpload(multiple=False, layout={'width': 'initial'})
        self.file_upload.observe(self._on_file_upload, names='value')
        self.callback_indicator = ipw.Valid(description='osp-callback')

        super().__init__(children=[self.file_upload, self.callback_indicator], **kwargs)
        
    @staticmethod
    def _to_ase(cuds_object):
        from osp.wrappers.simase import SimaseSession
        from osp.core.namespaces import simase_ontology

        session = SimaseSession()
        ase_wrapper = simase_ontology.simase_wrapper(session=session)
        ase_wrapper.add(cuds_object)
        session.run()
        return session._atoms
        
    def _on_file_upload(self, change=None):
        for fname, item in change['new'].items():
            assert item['metadata']['type'] == 'application/json'
            cuds_object = deserialize(json.loads(item['content'].decode('utf-8')))
            self.structure = self._to_ase(cuds_object)

    def from_callback(self, url):
        with self.hold_trait_notifications():
            try:
                query = parse_qs(urlsplit(url).query)
                osp_callback = query['osp_callback'][0]
            except (KeyError, IndexError):
                self.callback_indicator.value = False
            else:
                response = requests.get(osp_callback)
                response.raise_for_status()
                cuds_object = deserialize(response.json())
                self.structure = self._to_ase(cuds_object)
                self.callback_indicator.value = True

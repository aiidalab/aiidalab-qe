import json

import ase
import traitlets
import ipywidgets as ipw

from osp.core.utils import deserialize


class CUDSImporter(ipw.HBox):
    
    title = 'CUDS'
    structure = traitlets.Instance(ase.Atoms, allow_none=True)
    
    def __init__(self, **kwargs):
        self.file_upload = ipw.FileUpload(multiple=False, layout={'width': 'initial'})
        self.file_upload.observe(self._on_file_upload, names='value')
        super().__init__(children=[self.file_upload], **kwargs)
        
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

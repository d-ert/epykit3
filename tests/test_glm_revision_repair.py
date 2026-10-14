"""Changed GLM inference cannot reuse an old binary result fingerprint."""

from epykit import _glm
from epykit._dmc_config import DMCConfig


def test_glm_revision_changes_resume_fingerprint_and_records_contrast(monkeypatch):
    cfg = DMCConfig(test='glm')
    old = cfg.resume_signature_params(selected_test='glm', unite=False)
    lr = cfg.resume_signature_params(selected_test='lr', unite=False)
    monkeypatch.setattr(_glm, 'GLM_ALGORITHM_REVISION', 'test-next-revision', raising=False)
    new = cfg.resume_signature_params(selected_test='glm', unite=False)
    assert old != new
    assert lr == cfg.resume_signature_params(selected_test='lr', unite=False)
    for test_used, terms in [('glm', None), ('glm_contrast', ['Intercept', 'group'])]:
        meta = cfg.to_uns(test_used=test_used, n_sites=10, materialized=True,
                          unite=False, last_key='dmc', store_path=None, design_terms=terms)
        assert meta['glm_algorithm_revision'] == 'test-next-revision'

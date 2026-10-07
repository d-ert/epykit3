"""Public execution, provenance and refusal tests for the actual count engine."""
from pathlib import Path
import json
import numpy as np
import pytest
import epykit as ep
from epykit._dmc_engines import PUBLIC_ENGINES


def test_public_count_engine_runs_and_reports_its_model(synth_md_filtered):
    assert 'beta_binomial' in PUBLIC_ENGINES, 'The count engine is not a public choice'
    md = synth_md_filtered
    ep.tl.dmc(md,test='beta_binomial',beta_binomial_ci=False,tsv=False)
    frame = md.dmc
    assert frame is not None and len(frame)>0
    assert {'bb_rho','bb_lr','bb_df','log2_odds_ratio_bb'}.issubset(frame.columns)
    p = frame['pvalue'].to_numpy()
    assert np.all((p[np.isfinite(p)]>=0)&(p[np.isfinite(p)]<=1))
    assert frame['meth_diff_ci_lo'].is_nan().all()
    assert md.uns['dmc']['test_used'] == 'beta_binomial'
    assert md.uns['dmc']['statistical_model'] == 'per-replicate beta-binomial'
    assert md.uns['dmc']['beta_binomial_ci'] is False
    manifest = md.dmc_store._manifest
    entry = manifest['chroms'][0]
    assert 'prior_file' in entry and 'prior_sha256' in entry
    prior_path = md.dmc_store.path/entry['prior_file']
    assert prior_path.exists()
    prior = json.loads(prior_path.read_text())
    assert prior['log_rho_sd']>0
    assert prior['n_groups']>0
    assert prior['algorithm_revision'].startswith('count-likelihood-eb-bb-lognormal')


@pytest.mark.parametrize('knobs',[
    {'formula':'~ treatment','contrast':'treatment'},
    {'covariates':['age']}, {'smoothing':True}, {'use_smoothed':True},
    {'dispersion':'site'}, {'sep_fallback':True}, {'neighbour_combine':True},
    {'empirical_fdr':True},
])
def test_unsupported_requests_are_refused_before_running_an_old_engine(synth_md_filtered,knobs):
    assert 'beta_binomial' in PUBLIC_ENGINES, 'The count engine is not a public choice'
    md = synth_md_filtered
    with pytest.raises(ValueError,match='beta_binomial'):
        ep.tl.dmc(md,test='beta_binomial',tsv=False,**knobs)
    assert 'dmc' not in md.uns


def test_profile_interval_setting_invalidates_count_engine_cache(synth_md_filtered):
    assert 'beta_binomial' in PUBLIC_ENGINES, 'The count engine is not a public choice'
    md = synth_md_filtered
    ep.tl.dmc(md,test='beta_binomial',beta_binomial_ci=False,tsv=False)
    assert md.dmc['meth_diff_ci_lo'].is_nan().all()
    ep.tl.dmc(md,test='beta_binomial',beta_binomial_ci=True,tsv=False)
    f = md.dmc
    assert f['meth_diff_ci_lo'].is_finite().any()
    keep = f['pvalue'].is_finite().to_numpy()
    p = f['pvalue'].to_numpy()[keep]
    lo,hi = f['meth_diff_ci_lo'].to_numpy()[keep],f['meth_diff_ci_hi'].to_numpy()[keep]
    np.testing.assert_equal((lo>0)|(hi<0),p<.05)


def test_separate_output_stores_preserve_group_specific_prior_provenance(synth_md_filtered,tmp_path):
    from epykit.dmc import process_chromosomes_dmc
    import hashlib
    md = synth_md_filtered
    chrom = next(Path(md.store).glob('sample=*/chrom=*')).name.removeprefix('chrom=')
    case,control = list(md.treatment_ids),list(md.control_ids)
    def run(a,b,out):
        return process_chromosomes_dmc(md.store,a,b,test='beta_binomial',dispersion='eb',
            chromosomes=[chrom],beta_binomial_ci=False,out_dir=out,return_store=True)
    first = run(case,control,tmp_path/'first')
    entry = first._manifest['chroms'][0]
    assert 'prior_file' in entry, 'Learned prior must belong to its output store'
    prior_path = first.path/entry['prior_file']
    before = prior_path.read_bytes()
    case[0],control[0] = control[0],case[0]
    second = run(case,control,tmp_path/'second')
    assert prior_path.read_bytes() == before
    assert hashlib.sha256(before).hexdigest() == entry['prior_sha256']
    assert json.loads(before)['samples_case'] != json.loads(
        (second.path/second._manifest['chroms'][0]['prior_file']).read_text())['samples_case']
    assert not (Path(md.store)/'.bb_priors').exists()


def test_read_only_input_store_runs_with_a_writable_output(synth_md_filtered,tmp_path):
    from epykit.dmc import process_chromosomes_dmc
    import stat
    root = Path(synth_md_filtered.store)
    chrom = next(root.glob('sample=*/chrom=*')).name.removeprefix('chrom=')
    paths = [root]+list(root.rglob('*'))
    modes = {p:stat.S_IMODE(p.stat().st_mode) for p in paths}
    try:
        for p in paths:
            p.chmod(0o555 if p.is_dir() else 0o444)
        result = process_chromosomes_dmc(root,list(synth_md_filtered.treatment_ids),
            list(synth_md_filtered.control_ids),test='beta_binomial',dispersion='eb',
            chromosomes=[chrom],beta_binomial_ci=False,out_dir=tmp_path/'writable',return_store=True)
        assert len(result._manifest['chroms'])==1
        assert (result.path/result._manifest['chroms'][0]['prior_file']).is_file()
    finally:
        for p,mode in modes.items():
            p.chmod(mode)

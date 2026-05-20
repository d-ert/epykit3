# Distributed Backends

epykit supports distributed and GPU-accelerated execution for computationally
intensive analyses. The sequential (single-process) backend is the default and
is sufficient for most analyses. Distributed backends help when working with
whole-genome WGBS data at scale -- many samples, many chromosomes, or GLM
covariate models.

## Dask

Dask distributes chromosome-level tasks across a local or cluster-based worker
pool.

**Install:**

```bash
pip install 'epykit[distributed]'
```

**Usage:**

```python
import epykit as ep

md = ep.read_bismark("samples.csv", treatment_group="tumor",
                     control_group="normal", assembly="hg38")
ep.pp.filter_coverage(md)
ep.pp.normalize_coverage(md)
ep.pp.unite(md)

ep.tl.dmc(md, backend="dask", n_workers=4)
```

Each chromosome is submitted as an independent Dask task. With `n_workers=4`,
up to four chromosomes are processed in parallel. The Dask scheduler manages
memory and task scheduling automatically.

!!! tip "Cluster deployment"
    For cluster execution, start a Dask distributed client before calling
    `ep.tl.dmc()`. epykit will detect the existing client and submit tasks to
    it instead of creating a local cluster.

    ```python
    from dask.distributed import Client
    client = Client("scheduler-address:8786")
    ep.tl.dmc(md, backend="dask")
    ```

## Ray

Ray provides an alternative distributed backend with its own scheduler and
object store.

**Install:**

```bash
pip install 'epykit[ray]'
```

**Usage:**

```python
ep.tl.dmc(md, backend="ray", n_workers=4)
```

Like Dask, each chromosome is an independent Ray task. Ray's object store can
be advantageous when intermediate results are large, as it avoids redundant
serialization.

## GPU Acceleration (CuPy)

The GLM test backend (`test="glm"`) supports GPU-accelerated IRLS via CuPy.
This routes the iteratively reweighted least squares hot loop to the GPU.

**Install:**

```bash
pip install 'epykit[gpu]'
```

**Usage:**

```python
ep.tl.dmc(md, test="glm", formula="~ group + age", glm_backend="gpu")
```

The `glm_backend="gpu"` parameter only affects the `glm` test. Other tests
(`lr`, `score`, `fisher`, etc.) do not have a GPU path because their per-site
computation is not matrix-heavy enough to benefit from GPU offloading.

!!! note "CUDA requirement"
    CuPy requires a CUDA-capable GPU and a matching CUDA toolkit installation.
    See the [CuPy installation guide](https://docs.cupy.dev/en/stable/install.html)
    for details.

## GPU Acceleration (JAX)

JAX provides an alternative GPU backend for the GLM IRLS loop, with automatic
differentiation and XLA compilation.

**Install:**

```bash
pip install 'epykit[gpu_jax]'
```

**Usage:**

```python
ep.tl.dmc(md, test="glm", formula="~ group + age", glm_backend="gpu_jax")
```

JAX compiles the IRLS kernel with XLA on first invocation, which adds a
one-time compilation overhead. Subsequent chromosomes reuse the compiled kernel
and run faster than the CuPy path for large design matrices.

## When to Use Distributed Backends

| Scenario | Recommendation |
|----------|---------------|
| RRBS or targeted panels (< 1M sites) | Sequential is fast enough |
| WGBS with 2--6 samples, `lr` or `score` test | Sequential is usually fine (minutes) |
| WGBS with 10+ samples | Dask or Ray reduces wall time |
| GLM with multiple covariates on WGBS | Dask/Ray for chromosome parallelism, `gpu` for IRLS |
| Single chromosome or a few chromosomes | Sequential -- parallelism overhead exceeds benefit |
| Cluster or cloud environment | Dask with an external scheduler |

### Sequential Is the Default for a Reason

The sequential backend has zero overhead: no serialization, no inter-process
communication, no scheduler startup. For most analyses (RRBS, targeted panels,
WGBS with a handful of samples), it completes in minutes and is the simplest
to debug. Switch to a distributed backend only when wall time becomes a
bottleneck.

## Combining Distributed and GPU

The `backend` and `glm_backend` parameters are independent. You can combine
them:

```python
# Distribute chromosomes across 4 Dask workers, each using GPU for GLM IRLS
ep.tl.dmc(
    md,
    test="glm",
    formula="~ group + age + sex",
    backend="dask",
    n_workers=4,
    glm_backend="gpu",
)
```

This is the maximum-throughput configuration for large covariate-adjusted WGBS
analyses. Each Dask worker processes one chromosome at a time and offloads the
IRLS computation to the GPU.

## Installation Summary

| Extra | Command | What it provides |
|-------|---------|-----------------|
| `distributed` | `pip install 'epykit[distributed]'` | Dask backend |
| `ray` | `pip install 'epykit[ray]'` | Ray backend |
| `gpu` | `pip install 'epykit[gpu]'` | CuPy GLM acceleration |
| `gpu_jax` | `pip install 'epykit[gpu_jax]'` | JAX GLM acceleration |

Multiple extras can be installed together:

```bash
pip install 'epykit[distributed,gpu]'
```

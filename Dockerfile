FROM condaforge/mambaforge:24.9.2-0

WORKDIR /analysis
COPY environment.yml /analysis/environment.yml
RUN mamba env create -f /analysis/environment.yml && \
    mamba clean --all --yes

ENV CONDA_DEFAULT_ENV=lps-oantigen-analysis
ENV PATH=/opt/conda/envs/lps-oantigen-analysis/bin:$PATH
COPY . /analysis

CMD ["bash", "reproducibility/run_checks.sh"]

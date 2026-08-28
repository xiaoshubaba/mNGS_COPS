FROM mambaorg/micromamba:1.5.8

USER root
WORKDIR /opt/pops
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /opt/pops/environment.yml
COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml README.md USAGE.md POPS.main.py /opt/pops/
COPY --chown=$MAMBA_USER:$MAMBA_USER src /opt/pops/src
COPY --chown=$MAMBA_USER:$MAMBA_USER config /opt/pops/config
COPY --chown=$MAMBA_USER:$MAMBA_USER examples /opt/pops/examples

USER $MAMBA_USER
RUN micromamba create -y -n pops -f /opt/pops/environment.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/envs/pops/bin:$PATH
RUN pip install --no-deps /opt/pops

ENTRYPOINT ["pops"]
CMD ["--help"]

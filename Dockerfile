# NGS analysis environment
# Bioconda tools (samtools, bcftools, fastp, seqtk, ...) are Linux-only,
# so we build the conda env inside a Linux container.
FROM condaforge/miniforge3:latest

# Create the shared base environment from the versioned spec.
COPY shared/envs/base.yml /tmp/base.yml
RUN mamba env create -f /tmp/base.yml && mamba clean -a -y

# Make the env active for all subsequent commands and interactive shells.
ENV PATH=/opt/conda/envs/seq-base/bin:$PATH
RUN echo "conda activate seq-base" >> ~/.bashrc

WORKDIR /work
CMD ["bash"]

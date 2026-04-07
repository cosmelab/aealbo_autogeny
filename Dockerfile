# syntax=docker/dockerfile:1.4
# Albopictus Autogeny GWAS - Pixi-based Container
# Platform: linux/amd64 only
# Key: Set ENV PATH globally so singularity exec works

FROM --platform=linux/amd64 ghcr.io/prefix-dev/pixi:latest AS build

# CRITICAL: Set PATH globally for singularity exec compatibility
# Without this, python/R only available via pixi shell-hook (zsh only)
ENV PATH="/workspace/.pixi/envs/default/bin:/usr/local/bin:/usr/bin:$PATH"

# 1) Install system dependencies first
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    curl \
    wget \
    git \
    zsh \
    bash \
    unzip \
    # Graphics libraries for R plotting
    libpng-dev \
    libfreetype6-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libcairo2-dev \
    libtiff5-dev \
    libjpeg-dev \
    # Development libraries
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libbz2-dev \
    liblzma-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    # Java for snpEff
    default-jre \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# 2) Set up workspace and install dependencies via pixi
WORKDIR /workspace
COPY pixi.toml ./
RUN pixi install --frozen || pixi install

# 3) Install additional R packages via pixi task (OutFLANK, LDna from GitHub, etc.)
RUN pixi run install-r-extras || echo "Some R packages may need manual installation"

# 4) Install Oh My Zsh and shell environment
RUN git clone --depth=1 https://github.com/ohmyzsh/ohmyzsh.git ~/.oh-my-zsh && \
    git clone --depth=1 https://github.com/romkatv/powerlevel10k.git ~/.oh-my-zsh/custom/themes/powerlevel10k && \
    git clone --depth=1 https://github.com/zsh-users/zsh-completions.git ~/.oh-my-zsh/custom/plugins/zsh-completions && \
    git clone --depth=1 https://github.com/zsh-users/zsh-autosuggestions.git ~/.oh-my-zsh/custom/plugins/zsh-autosuggestions && \
    git clone --depth=1 https://github.com/zsh-users/zsh-syntax-highlighting.git ~/.oh-my-zsh/custom/plugins/zsh-syntax-highlighting && \
    cp ~/.oh-my-zsh/templates/zshrc.zsh-template ~/.zshrc && \
    sed -i 's/ZSH_THEME=".*"/ZSH_THEME="powerlevel10k\/powerlevel10k"/' ~/.zshrc && \
    sed -i 's/plugins=(git)/plugins=(git zsh-completions zsh-autosuggestions zsh-syntax-highlighting)/' ~/.zshrc && \
    echo 'export DISABLE_AUTO_UPDATE="true"' >> ~/.zshrc && \
    echo 'export DISABLE_UPDATE_PROMPT="true"' >> ~/.zshrc && \
    echo 'POWERLEVEL9K_DISABLE_CONFIGURATION_WIZARD=true' >> ~/.zshrc

# 5) Install fzf
RUN mkdir -p ~/.fzf && \
    git clone --depth 1 https://github.com/junegunn/fzf.git ~/.fzf && \
    ~/.fzf/install --all && \
    echo '[ -f ~/.fzf.zsh ] && source ~/.fzf.zsh' >> ~/.zshrc && \
    echo 'export FZF_BASE=~/.fzf' >> ~/.zshrc

# 6) Configure shell environment
RUN echo 'if command -v eza > /dev/null; then' >> ~/.zshrc && \
    echo '  alias ls="eza"' >> ~/.zshrc && \
    echo '  alias ll="eza -l"' >> ~/.zshrc && \
    echo '  alias la="eza -la"' >> ~/.zshrc && \
    echo '  alias lt="eza --tree"' >> ~/.zshrc && \
    echo 'fi' >> ~/.zshrc && \
    echo 'export TERM="xterm-256color"' >> ~/.zshrc && \
    echo 'cd /workspace && eval "$(pixi shell-hook)"' >> ~/.zshrc && \
    echo 'eval "$(starship init zsh)"' >> ~/.zshrc

# 7) Configure bash fallback
RUN echo 'cd /workspace && eval "$(pixi shell-hook)"' >> ~/.bashrc

# 8) Set up workspace directory structure
WORKDIR /workspace
RUN mkdir -p /workspace/data/{raw,metadata,references,genome,genotype_calls,files} \
             /workspace/output/{quality_control,outflank,pcadapt,ldna,snpeff,fst,missingness} \
             /workspace/scripts \
             /workspace/notebooks \
             /workspace/logs \
             /workspace/results/{figures,tables,supplementary}

# 9) Copy project files
COPY . /workspace

# 10) Make scripts executable and create R profile
RUN find /workspace/scripts -name "*.py" -exec chmod +x {} \; 2>/dev/null || true && \
    find /workspace/scripts -name "*.sh" -exec chmod +x {} \; 2>/dev/null || true && \
    echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"))' > /workspace/.Rprofile && \
    echo 'options(download.file.method = "libcurl")' >> /workspace/.Rprofile && \
    echo 'options(timeout = 600)' >> /workspace/.Rprofile

# 11) Add welcome message
RUN echo '' >> ~/.zshrc && \
    echo 'echo "==================================================="' >> ~/.zshrc && \
    echo 'echo "==================================================="' >> ~/.zshrc && \
    echo 'echo "   Albopictus Autogeny Analysis"' >> ~/.zshrc && \
    echo 'echo "   Reproducible R notebooks + CLI pipeline"' >> ~/.zshrc && \
    echo 'echo "==================================================="' >> ~/.zshrc && \
    echo 'echo "Genomics: plink, plink2, bcftools, vcftools, snpEff"' >> ~/.zshrc && \
    echo 'echo "R: tidyverse, OutFLANK, pcadapt, adegenet, vcfR, LDna"' >> ~/.zshrc && \
    echo 'echo "Render:  Rscript -e rmarkdown::render(notebook.Rmd)"' >> ~/.zshrc && \
    echo 'echo ""' >> ~/.zshrc && \
    echo 'echo "Run: pixi run test-tools  # verify installation"' >> ~/.zshrc && \
    echo 'echo "==================================================="' >> ~/.zshrc

# 12) Verify key tools
RUN echo "Verifying installations..." && \
    pixi run test-tools || echo "Some tools may need verification"

ENV R_PROFILE_USER=/workspace/.Rprofile
WORKDIR /workspace

ENTRYPOINT ["zsh", "-l", "-c", "source ~/.zshrc && exec zsh"]
CMD []

# NGS workshop login environment, installed as /etc/profile.d/ngs.sh by
# build_image.sh. @REF@ and @GFF3@ are filled in at build time.
if [ -n "${BASH_VERSION:-}" ] && [ -f /opt/miniforge/etc/profile.d/conda.sh ]; then
  . /opt/miniforge/etc/profile.d/conda.sh
  conda activate ngs
fi
export PATH="/opt/ngs/bin:$PATH"
export NGS_REF=@REF@
export FASTVEP_GFF3=@GFF3@
export FASTVEP_SA_DIR=/opt/assets/fastvep/sa

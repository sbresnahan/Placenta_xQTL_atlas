# Porting to another system

- Replace seadragon paths (`/rsrch5`, `/rsrch9`) and infrastructure (conda
  env names/paths, singularity image, `R_LIBS_USER`, module loads). These
  are hardcoded in scripts; data and reference paths are config-driven via
  `03_phenotyping/config.yml`.
- LSF to SLURM: convert `#BSUB` to `#SBATCH`, `bsub` to `sbatch`,
  `-w "done(...)"` to `--dependency=afterok:...`, and the submission calls
  in `03_phenotyping/run_pipeline.py`.
- Any R >= 4.3 with the packages in
  [software_environments.md](software_environments.md) works; the singularity
  image is only needed to reproduce the exact seadragon stack.

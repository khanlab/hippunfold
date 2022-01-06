This workflow generates a set of custom density unfolded templates, with adaptive spacing designed to result in (mostly) constant surface area when transformed to native space



Requirements:
  - https://github.com/gllmflndn/gifti (clone into this folder)
  - matlab
  - hippunfold with uniform template run on a set of subjects
    - here we used HCP UR100, and hippunfold unfoldDG branch
    - it uses this to create an average surface area (gyrification) map

To run for hipp surf:
```
snakemake -np --configfile config_dentate.yml
```

To run for dentate surf:
```
snakemake -np --configfile config_dentate.yml
```
